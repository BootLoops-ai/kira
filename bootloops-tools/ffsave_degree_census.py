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
r"""
ffsave_degree_census.py — decode a (possibly still-running / never-finishing) kira+FireFly
reduction's ff_save/ state into a per-entry (target, master) -> (max_deg_num, max_deg_den,
is_done) census, WITHOUT any new solver compute.

Why: FireFly's first-prime degree scan measures every requested entry's rational-function
degree long before reconstruction finishes. On runs that will never finish (deg~300+
entries needing many primes), the degrees are the *diagnostic deliverable* — they tell you
WHICH masters/targets carry a degree explosion (basis-rotation targeting, target bisection)
for free. Validated on a long-running reduction (6793 entries):
decode gate = every tag component must land in the computed weight set of the KNOWN
de_targets / masters.final names (exact set match, 6793/6793).

Kira weight encoding (derived from pyred/integrals.cpp and verified empirically —
unique layout survived a 768-layout brute force; codec shared via pyred_weight.py):
  weight = (((((topo << np | sector_w) << b_d1 | d1) << b_d2 | d2)
             << b_sp | spows_w) << b_pp | ppows_w)
  - sector_w: sector rank. integral_ordering 1-4 => sector_ordering 1 (rank =
    sector number); 5-8 => sector_ordering 2 (rank in sort by (#lines, sector
    number)) — the mapping is kira/ReadYamlFiles.cpp's io -> (so, dotsp).
  - (d1, d2) per dotsp ordering (= io for 1-4, io-4 for 5-8):
      1: (dots, sps)   2: (sps, dots)   3: (dots+sps, dots)   4: (dots+sps, sps)
    so integral_ordering 8 has d1 = dots+sps and integral_ordering 5 has d1 = dots.
    (ordering-5 layout: unique survivor of a 32-layout brute force gated on three
     independent 4-loop tag sets [122+234+64 targets, 457+583+476 masters, exact];
     cross-validated against independently derived weights: a pinned target weight
     reproduced exactly, 20/20 sector-rank spots, 1356/1356 exact name->weight
     matches on singleton EXACT-SET classes. Orderings 5 and 8 are
     validated on real reductions; the remaining branches are transcribed from
     pyred/integrals.cpp weight_to_dots/sector_weight_table and guarded by the
     same round-trip layout gate.)
  - ppows_w = index of the dot composition in kira's compositions(dots, lines) enumeration
  - spows_w = (#compositions(sps, nums) - 1) - index of the sp composition (reversed)
  - bit widths (b_d1, b_d2, b_pp, b_sp) = WEIGHTBITS row of results/kira.db
  - preferred_masters get custom weights 1..N in FILE ORDER, skipping trivial-sector entries
Tags in ff_save/tags (and states/*.gz tag_name) are "<targetWeight>_<masterWeight>".

Multivariate runs (kira run_firefly with symbolic invariants): states carry
individual_degrees_num/den, one per FireFly-internal variable position; ff_save/var_order.gz
maps position -> kira variable index (kira hands FireFly [d, <kinematic_invariants in
kinematics.yaml order>] as x1..xn; FireFly may permute — the log line "Using optimized
variable order" and var_order.gz record the permutation). Pass --vars d,s2,s12 (the kira
x1..xn names) to get named per-variable degrees; mapping is recorded in the output either
way so entries can be relabeled without recompute.

SAVE-FORMAT VARIANT (seen on eta-chain runs): some runs write the
individual_degrees_num/den sections EMPTY and no ff_save/var_order.gz ("Using default
variable order" in firefly.log => identity permutation). Per-variable degrees are still
fully determined by the monomial EXPONENT VECTORS stored in the g_ni/g_di sections
(reconstructed rational coefficients, lines "e1..ek <num> <den>") and combined_ni/
combined_di (unstable CRT residues, lines "e1..ek <residue>"), as FireFly's
RatReconst::save_state writes them (decoder validated against an independently harvested per-variable
census, 27,068 fns exact). With --vars, empty individual_degrees now
fall back to that harvest automatically; each entry carries ind_src
("individual_degrees" | "monomial_exponents" | "none") so consumers know the source.
NOTE: the fallback needs k = len(--vars) to parse exponent columns — it is only active
with --vars. ZERO-state files (literal first line "ZERO") are counted and skipped.

Second gap in the same class: eta-chain/DE-seed runs tag integrals that
live in NO run file (dotted masters, corner seeds) — name-list matching cannot cover
them, and the trivial-sector auto-search is exponential in distinct preferred sectors
(2^181 on one run). Cures: (1) the weight encoding is a bijection, so entries are now
named by INVERSE decode (weight -> name directly) whenever the file maps miss; per-entry
target_src/master_src say which path named it ('file' | 'custom' | 'inverse'). (2) the
layout gate is the round-trip: every file name must satisfy inverse(encode(name)) ==
name (assert), and every observed weight must decode (fail-loud ValueError) — same
catches-misuse philosophy, no exponential search. (3) the auto-search is capped at
r <= 2 dropped sectors; when no small hypothesis closes the old file gate the run
reports mode=inverse-decode instead of dying (file_gate field in the output).

usage: ffsave_degree_census.py <rundir> [--family FAM] [--np 0] [--vars d,s2,s12]
                               [--workers 1] [--out census.json]
       <rundir> must contain ff_save/states/, de_targets (or target), preferred,
       results/<family>/masters, results/kira.db (for WEIGHTBITS + INTEGRALORDERING).
       --np 0 (default) infers the propagator count from the name tuples.
NOTE: snapshot ff_save/states first if kira is live (files rotate).
Implements kira's full integral_ordering range 1..8 (see pyred_weight.py for the
per-ordering layout provenance); anything outside that range fails closed — the
round-trip layout gate catches a wrong-layout misread either way.
"""
import os, sys, gzip, re, json, sqlite3, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pyred_weight import WeightCodec

ap = argparse.ArgumentParser()
ap.add_argument('rundir')
ap.add_argument('--family', default='gtw')
ap.add_argument('--np', type=int, default=0,
                help='propagator count; 0 = infer from name tuple length')
ap.add_argument('--vars', default=None,
                help='comma list naming kira FireFly vars x1..xn (e.g. d,s2,s12); '
                     'enables named per-variable degrees via ff_save/var_order.gz')
ap.add_argument('--out', default=None)
ap.add_argument('--workers', type=int, default=1,
                help='parallel state-file readers (gzip-bound; 16 turns a ~45-60 min '
                     '599M/27k-file census into minutes). Default 1 = serial.')
ap.add_argument('--trivial-sectors', default='auto',
                help='comma list, or "auto": infer by dropping preferred entries whose '
                     'exclusion makes the tag gate pass')
a = ap.parse_args()
D, FAM = a.rundir, a.family

con = sqlite3.connect(f'file:{os.path.join(D,"results","kira.db")}?mode=ro&immutable=1', uri=True)
(b_d1, b_d2, b_pp, b_sp), = con.execute('SELECT * FROM WEIGHTBITS').fetchall()
(iord,), = con.execute('SELECT * FROM INTEGRALORDERING').fetchall()
con.close()
assert 1 <= iord <= 8, \
    f'integral_ordering {iord}: engine range is 1..8'

def names(fp):
    out = []
    for ln in open(fp):
        m = re.match(re.escape(FAM) + r'\[([^\]]+)\]', ln.strip())
        if m: out.append(tuple(int(x) for x in m.group(1).split(',')))
    return out

masters = names(os.path.join(D, 'results', FAM, 'masters'))
preferred = names(os.path.join(D, 'preferred'))
tfile = os.path.join(D, 'de_targets')
if not os.path.exists(tfile): tfile = os.path.join(D, 'target')
targets = names(tfile)

NP = a.np or max(len(nu) for nu in masters + targets)
assert all(len(nu) == NP for nu in masters + targets), 'mixed name tuple lengths; pass --np'

CODEC = WeightCodec(NP, (b_d1, b_d2, b_pp, b_sp), iord)
CUSTOM_CEIL = 1 << (b_d1 + b_d2 + b_sp + b_pp)   # any default weight of a sector
# rank >= 1 exceeds this; observed weights below it are unambiguously custom
def sec_of(nu): return sum(1 << i for i, p in enumerate(nu) if p > 0)
def default_weight(nu):
    """Name -> default pyred weight, per the run's integral_ordering."""
    return CODEC.encode(nu)

def inverse_weight(w):
    """Weight -> integral name, the layout run backwards (DE-seed tags name
    dotted/corner integrals that live in NO run file, so name-list matching
    cannot cover them — but the
    weight encoding is a bijection, so decode directly). Raises ValueError with
    the failing component on any out-of-range index (fail loud, never guess)."""
    return CODEC.decode(w)

KEYS = ('tag_name', 'is_done', 'max_deg_num', 'max_deg_den',
        'individual_degrees_num', 'individual_degrees_den')
ALLKEYS = KEYS + ('combined_prime', 'need_prime_shift', 'normalizer_deg',
                  'normalize_to_den', 'normalizer_den_num',
                  'shifted_max_num_eqn', 'shift', 'sub_num', 'sub_den',
                  'zero_degs_num', 'zero_degs_den', 'ni', 'di')
# every section keyword FireFly's save_state writes — needed to
# delimit the monomial sections during the exponent harvest
SECKEYS = frozenset(ALLKEYS + ('g_ni', 'g_di', 'combined_ni', 'combined_di',
                               'combined_primes_ni', 'combined_primes_di',
                               'interpolations'))
NVARS = len(a.vars.split(',')) if a.vars else 0

def pervar_from_monomials(lines):
    """Harvest per-variable degree maxima from the monomial exponent vectors of
    g_ni/g_di (fields: e1..ek num den) and combined_ni/combined_di (fields:
    e1..ek residue). Returns (ind_n, ind_d) or None if no monomial line found.
    Field counts are asserted — a mismatch means wrong --vars arity or an
    unimplemented format, and the gate philosophy is fail-loud, never guess."""
    MSECS = {'g_ni': (0, NVARS + 2), 'combined_ni': (0, NVARS + 1),
             'g_di': (1, NVARS + 2), 'combined_di': (1, NVARS + 1)}
    mx = [[0] * NVARS, [0] * NVARS]; seen = False
    cur = None
    for l in lines:
        if l in SECKEYS:
            cur = l if l in MSECS else None; continue
        if not l or cur is None: continue
        side, nf = MSECS[cur]
        parts = l.split()
        assert len(parts) == nf, (
            f'monomial line in {cur} has {len(parts)} fields, expected {nf} '
            f'for {NVARS} vars — wrong --vars arity or unimplemented format')
        seen = True
        for v in range(NVARS):
            e = int(parts[v])
            if e > mx[side][v]: mx[side][v] = e
    return (mx[0], mx[1]) if seen else None

def parse_state(fn):
    """One states/*.gz -> ((w1, w2), rec) | ('ZERO', fn) | None (unreadable)."""
    try:
        with gzip.open(os.path.join(sdir, fn), 'rt') as f: lines = f.read().split('\n')
    except (FileNotFoundError, EOFError, OSError): return None
    if lines and lines[0] == 'ZERO':
        return ('ZERO', fn)           # zero function (save_zero_state); counted, no degrees
    kv = {}
    for i, l in enumerate(lines[:60]):
        if l in KEYS:
            v = lines[i + 1] if i + 1 < len(lines) else ''
            # univariate runs: a key with an EMPTY value is immediately followed by
            # the next key line — do not eat the key name as a value
            kv[l] = '' if v.strip() in ALLKEYS else v
    w1, w2 = (int(x) for x in kv['tag_name'].split('_'))
    ind_n = [int(x) for x in kv.get('individual_degrees_num', '').split()]
    ind_d = [int(x) for x in kv.get('individual_degrees_den', '').split()]
    ind_src = 'individual_degrees' if (ind_n or ind_d) else 'none'
    if a.vars and not ind_n and not ind_d:
        # save-format variant: empty individual_degrees sections; the
        # per-var maxima live in the monomial exponent vectors instead
        got = pervar_from_monomials(lines)
        if got is not None:
            ind_n, ind_d = got; ind_src = 'monomial_exponents'
    # short-format states (post-degree-scan, pre-interpolation) carry only tag_name,
    # normalize_to_den and the individual degrees: max_deg_* then falls back to the
    # per-variable maxima (a lower bound on the total degree; flagged via has_max=0).
    has_max = int('max_deg_num' in kv)
    dn = int(kv['max_deg_num']) if has_max else (max(ind_n) if ind_n else 0)
    dd = int(kv['max_deg_den']) if has_max else (max(ind_d) if ind_d else 0)
    return ((w1, w2), (dn, dd, int(kv.get('is_done', 0)), ind_n, ind_d, has_max, ind_src))

entries = {}; nzero = 0
sdir = os.path.join(D, 'ff_save', 'states')
_files = os.listdir(sdir)
if a.workers > 1:
    import multiprocessing as _mp
    if _mp.get_start_method() != 'fork':
        print(f"[census] WARN: --workers>1 needs the fork start method "
              f"(Linux default); '{_mp.get_start_method()}' (macOS default) "
              f"re-runs this script's top level in every worker — "
              f"running SERIAL instead", file=sys.stderr)
        a.workers = 1
if a.workers > 1:
    from multiprocessing import Pool
    with Pool(a.workers) as _pool:
        _parsed = _pool.map(parse_state, _files, chunksize=64)
else:
    _parsed = map(parse_state, _files)
for r in _parsed:
    if r is None: continue
    if r[0] == 'ZERO': nzero += 1; continue
    key, rec = r
    old = entries.get(key)
    if old is None or (rec[2], rec[0] + rec[1]) > (old[2], old[0] + old[1]):
        entries[key] = rec            # keep most-advanced state if a tag has several files

# per-variable position -> name mapping (multivariate runs)
var_names = None; var_order = None
vo_fp = os.path.join(D, 'ff_save', 'var_order.gz')
if os.path.exists(vo_fp):
    with gzip.open(vo_fp, 'rt') as f:
        var_order = {int(p): int(v) for p, v in
                     (ln.split() for ln in f.read().strip().split('\n') if ln.strip())}
if a.vars:
    kv_names = a.vars.split(',')
    if var_order:
        var_names = [kv_names[var_order[p]] for p in sorted(var_order)]
    else:
        var_names = kv_names

# G2 — LAYOUT GATE: every file name must round-trip encode->decode exactly.
# This is the check that catches a wrong ordering/layout (the old exact-set
# gate's real content), independent of whether the run's tags stay inside the
# file name lists.
for nu in set(masters) | set(targets) | set(preferred):
    back = inverse_weight(default_weight(nu))
    assert back == nu, f'LAYOUT GATE FAIL: {nu} -> {default_weight(nu)} -> {back}'

def build_maps(trivial):
    custom = {}; k = 0
    for nu in preferred:
        if sec_of(nu) in trivial: continue
        k += 1; custom[nu] = k
    wfun = lambda nu: custom.get(nu) or default_weight(nu)
    return {wfun(m): m for m in masters}, {wfun(t): t for t in targets}

W1 = {k[0] for k in entries}; W2 = {k[1] for k in entries}
AUTO_MAX_R = 2   # an unbounded combinations search is exponential in the
# number of distinct preferred sectors (2^181 on one production run — a wall,
# not a wait); realistic trivial sets are tiny, so cap it
# and fall back to inverse decode when no small hypothesis closes the file gate.
trivial = set(); file_gate = False
if a.trivial_sectors == 'auto':
    from itertools import combinations
    psecs = sorted({sec_of(nu) for nu in preferred})
    for r in range(min(len(psecs), AUTO_MAX_R) + 1):
        hit = next((set(cmb) for cmb in combinations(psecs, r)
                    if (lambda mw_tw: W2 <= set(mw_tw[0]) and W1 <= set(mw_tw[1]))(build_maps(set(cmb)))), None)
        if hit is not None: trivial = hit; file_gate = True; break
else:
    trivial = {int(x) for x in a.trivial_sectors.split(',')} if a.trivial_sectors else set()
    mw0, tw0 = build_maps(trivial)
    file_gate = W2 <= set(mw0) and W1 <= set(tw0)
mw, tw = build_maps(trivial)

# name resolution: file maps first (old behavior, byte-compatible when the file
# gate closes), custom weights via surviving-preferred order, else INVERSE decode
customs_by_w = {}; k = 0
for nu in preferred:
    if sec_of(nu) in trivial: continue
    k += 1; customs_by_w[k] = nu
def resolve(w, filemap):
    if w in filemap: return filemap[w], 'file'
    if w < CUSTOM_CEIL:
        nu = customs_by_w.get(w)
        if nu is None:
            raise ValueError(f'custom weight {w} exceeds the surviving preferred list '
                             f'({len(customs_by_w)}) — trivial-sector hypothesis wrong?')
        return nu, 'custom'
    return inverse_weight(w), 'inverse'   # G1: raises on any undecodable weight

resolved = {(w1, w2): (resolve(w1, tw), resolve(w2, mw)) for (w1, w2) in entries}
nsrc = {}
for (t_, m_) in resolved.values():
    for _, s in (t_, m_): nsrc[s] = nsrc.get(s, 0) + 1
srcs = {}
for rec in entries.values(): srcs[rec[6]] = srcs.get(rec[6], 0) + 1
mode = ('file-complete' if file_gate else
        'inverse-decode (tags name integrals outside the run name files — DE-seed class)')
print(f'[ffdeg] GATE PASS [{mode}]: {len(entries)} entries decode; trivial sectors '
      f'{sorted(trivial)}; {len(W1)} targets x {len(W2)} masters; name sources {nsrc}; '
      f'zero-states {nzero}; pervar sources {srcs}', file=sys.stderr)
out = {'rundir': D, 'weight_bits': [b_d1, b_d2, b_pp, b_sp], 'integral_ordering': iord,
       'np': NP, 'trivial_sectors': sorted(trivial), 'file_gate': file_gate,
       'var_order': var_order, 'var_names': var_names,
       'zero_states': nzero, 'pervar_sources': srcs, 'name_sources': nsrc,
       'entries': [{'tag': f'{w1}_{w2}',
                    'target': list(resolved[(w1, w2)][0][0]),
                    'master': list(resolved[(w1, w2)][1][0]),
                    'target_src': resolved[(w1, w2)][0][1],
                    'master_src': resolved[(w1, w2)][1][1],
                    'deg_num': dn, 'deg_den': dd, 'done': done,
                    'ind_num': ind_n, 'ind_den': ind_d, 'has_max': hm,
                    'ind_src': isrc}
                   for (w1, w2), (dn, dd, done, ind_n, ind_d, hm, isrc) in entries.items()]}
json.dump(out, open(a.out or os.path.join(D, 'FF_DEGREE_CENSUS.json'), 'w'))
print(f'[ffdeg] wrote {a.out or os.path.join(D, "FF_DEGREE_CENSUS.json")}', file=sys.stderr)
