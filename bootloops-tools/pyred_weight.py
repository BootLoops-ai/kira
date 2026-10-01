#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Derived from Kira (https://gitlab.com/kira-pyred/kira), Copyright (C)
# 2017-2025 The Kira Developers (see ../AUTHORS): the integral-weight layout
# follows pyred/integrals.cpp (Topology::sector_weight_table,
# Integral::to_weight, Integral::weight_to_dots) and the integral_ordering
# mapping of kira/ReadYamlFiles.cpp; the ff_save reader follows FireFly's
# on-disk state format (FireFly, GPL-3.0, J. Klappert, S. Y. Klein, F. Lange).
# Python port and extensions Copyright (c) 2026 Anthropic, PBC; created by
# Matthew D. Schwartz, code written by Claude (Anthropic) under his supervision.
# This file is distributed under the GNU General Public License v3.0 or
# later, the license of the works it derives from; see ../COPYING.
r"""pyred_weight — Kira/pyred integral-weight codec, integral_ordering 1..8.

Encodes integral names (index tuples) to pyred weight IDs and decodes weights
back to names, matching the engine layout exactly. The semantics are
transcribed from the engine sources in this repository (src/kira, src/pyred):

  integral_ordering -> (sector_ordering, dotsp_ordering)
      kira/ReadYamlFiles.cpp: orderings 1..4 use sector_ordering 1 with
      dotsp_ordering = io; orderings 5..8 use sector_ordering 2 with
      dotsp_ordering = io - 4.
  sector rank
      pyred/integrals.cpp Topology::sector_weight_table:
      sector_ordering 1 ranks by sector number alone; sector_ordering 2
      ranks by (#lines, sector number), with sectors that are NOT a
      subsector of any declared top sector ordered after in-tree sectors
      of the same line count (the 1<<np tie-break offset).
  dotsp fields and composition indices
      pyred/integrals.cpp Integral::to_weight / weight_to_dots.

Layout (top bits to bottom; the WEIGHTBITS row of results/kira.db gives
b_d1, b_d2, b_pp, b_sp):

  weight = ((((rank << b_d1 | d1) << b_d2 | d2) << b_sp | spw) << b_pp | ppw)

with (d1, d2) per dotsp_ordering:
  1: (dots, sps)    2: (sps, dots)    3: (dots+sps, dots)    4: (dots+sps, sps)

ppw = index of the dot composition in the compositions(dots, lines)
enumeration; spw = (#compositions(sps, nums) - 1) - index of the sp
composition (reversed) — sp compositions to the left are lower weight,
dot compositions to the right are lower weight, exactly as in to_weight.

Assumptions shared with the rest of this package: topology id 0
(single-family run dirs) and an identity propagator permutation.
Preferred-master custom weights (small ints in file order) are a separate
name space layered on top by the callers — see ffsave_degree_census.py here
and the top-emit route of parallel_kira_gen.py in the bootloops toolkit
(tools/kira-stack).

Fail-loud: decoding raises ValueError with the failing component on any
out-of-range index; it never guesses. Encoding asserts every field fits its
declared bit width.
"""
import sqlite3
from functools import lru_cache


@lru_cache(maxsize=None)
def comps(n, k):
    """All weak compositions of n into k parts, in pyred enumeration order."""
    if k == 0:
        return ((),) if n == 0 else ()
    compos = [[]]
    for i in range(1, k):
        tmp = []
        for c in compos:
            cl = sum(c)
            for j in range(n + 1 - cl):
                tmp.append(c + [j])
        compos = tmp
    return tuple(tuple(c + [n - sum(c)]) for c in compos)


@lru_cache(maxsize=None)
def cidx(n, k):
    return {c: i for i, c in enumerate(comps(n, k))}


def sec_of(nu):
    """Sector bitmask of a name tuple (bit i set iff index i is positive)."""
    return sum(1 << i for i, p in enumerate(nu) if p > 0)


def from_kiradb(dbpath):
    """results/kira.db -> ((b_d1, b_d2, b_pp, b_sp), integral_ordering)."""
    con = sqlite3.connect(f'file:{dbpath}?mode=ro&immutable=1', uri=True)
    (bits,) = con.execute('SELECT * FROM WEIGHTBITS').fetchall()
    ((iord,),) = con.execute('SELECT * FROM INTEGRALORDERING').fetchall()
    con.close()
    return tuple(bits), iord


class WeightCodec:
    """Encode/decode pyred default weights for one run configuration.

    np:         number of indices (propagators + ISPs) in the family.
    bits:       (b_d1, b_d2, b_pp, b_sp) — the WEIGHTBITS row.
    iord:       kira integral_ordering, 1..8.
    topsectors: declared top sectors of the family (sector_ordering 2 ranks
                out-of-tree sectors after in-tree ones of the same line
                count). Default: the all-lines sector — every sector is
                in-tree, which is also the behavior validated on production
                censuses.
    """

    def __init__(self, np, bits, iord, topsectors=None):
        if not 1 <= iord <= 8:
            raise ValueError(
                f'integral_ordering {iord}: engine range is 1..8')
        if np < 2:
            raise ValueError('np=1 uses a special-cased layout '
                             '(Integral::to_weight) — not implemented')
        self.np = np
        self.b_d1, self.b_d2, self.b_pp, self.b_sp = bits
        self.iord = iord
        self.sector_ordering = 1 if iord <= 4 else 2
        self.dotsp = iord if iord <= 4 else iord - 4
        tops = list(topsectors) if topsectors else [(1 << np) - 1]
        if self.sector_ordering == 1:
            self.sects = list(range(1 << np))
        else:
            def key(s):
                intree = any((s & t) == s for t in tops)
                return (bin(s).count('1'), s + (0 if intree else 1 << np))
            self.sects = sorted(range(1 << np), key=key)
        self.rank = {s: i for i, s in enumerate(self.sects)}

    # ------------------------------------------------------------- encode

    def _enc_dotsp(self, dots, sps):
        o = self.dotsp
        if o == 1: return dots, sps
        if o == 2: return sps, dots
        if o == 3: return dots + sps, dots
        return dots + sps, sps

    def encode(self, nu):
        """Name tuple -> default pyred weight (topology id 0)."""
        assert len(nu) == self.np, f'{nu}: expected {self.np} indices'
        dots = sps = 0
        pp = []
        sp = []
        for p in nu:
            if p > 0:
                dots += p - 1
                pp.append(p - 1)
            else:
                sps -= p
                sp.append(-p)
        s = sec_of(nu)
        lines = bin(s).count('1')
        nums = self.np - lines
        ppw = cidx(dots, lines)[tuple(pp)]
        spw = len(comps(sps, nums)) - 1 - cidx(sps, nums)[tuple(sp)]
        d1, d2 = self._enc_dotsp(dots, sps)
        for v, b, what in ((d1, self.b_d1, 'd1'), (d2, self.b_d2, 'd2'),
                           (spw, self.b_sp, 'spw'), (ppw, self.b_pp, 'ppw')):
            assert v < (1 << b), f'{nu}: {what}={v} exceeds {b} weight bits'
        w = self.rank[s]
        w = (w << self.b_d1) | d1
        w = (w << self.b_d2) | d2
        w = (w << self.b_sp) | spw
        w = (w << self.b_pp) | ppw
        return w

    # ------------------------------------------------------------- decode

    def _dec_dotsp(self, w, d1, d2):
        o = self.dotsp
        if o == 1: dots, sps = d1, d2
        elif o == 2: sps, dots = d1, d2
        elif o == 3: dots, sps = d2, d1 - d2
        else: dots, sps = d1 - d2, d2
        if dots < 0 or sps < 0:
            raise ValueError(f'weight {w}: negative dots/sps '
                             f'({d1=} {d2=} dotsp_ordering {o})')
        return dots, sps

    def decode(self, w):
        """Weight -> name tuple, the layout run backwards (bijective on
        default weights). Raises ValueError with the failing component on
        any out-of-range index — fail loud, never guess."""
        ppw = w & ((1 << self.b_pp) - 1); r = w >> self.b_pp
        spw = r & ((1 << self.b_sp) - 1); r >>= self.b_sp
        d2 = r & ((1 << self.b_d2) - 1); r >>= self.b_d2
        d1 = r & ((1 << self.b_d1) - 1); rank = r >> self.b_d1
        if not (0 <= rank < len(self.sects)):
            raise ValueError(f'weight {w}: sector rank {rank} out of range')
        s = self.sects[rank]
        lines = bin(s).count('1')
        nums = self.np - lines
        dots, sps = self._dec_dotsp(w, d1, d2)
        cpp = comps(dots, lines)
        csp = comps(sps, nums)
        if not (0 <= ppw < len(cpp)):
            raise ValueError(f'weight {w}: ppw {ppw} out of range '
                             f'({len(cpp)} comps)')
        isp = len(csp) - 1 - spw
        if not (0 <= isp < len(csp)):
            raise ValueError(f'weight {w}: spw {spw} out of range '
                             f'({len(csp)} comps)')
        pp = cpp[ppw]
        sp = csp[isp]
        nu = []
        ip = isn = 0
        for i in range(self.np):
            if (s >> i) & 1:
                nu.append(1 + pp[ip]); ip += 1
            else:
                nu.append(-sp[isn]); isn += 1
        return tuple(nu)
