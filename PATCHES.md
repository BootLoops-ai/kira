# Kira fork — parallel-treatcoeff + SQLite autocommit guards (with kira128 / jemalloc build variants)

This repository is a lightly patched fork of **Kira**, the integration-by-parts
(IBP) reduction program (Laporta's algorithm with finite-field and algebraic
back ends): a snapshot at upstream release 3.1 with the two BootLoops patches
below applied; `PATCHES.diff` reproduces the complete change relative to that
base. Upstream's history and releases live at gitlab.com/kira-pyred/kira.

- **Upstream:** [gitlab.com/kira-pyred/kira](https://gitlab.com/kira-pyred/kira),
  by the Kira developers — created by Philipp Maierhöfer, Johann Usovitsch and
  Peter Uwer (Kira 1.0), joined by Jonas Klappert and Fabian Lange (Kira 2.0)
  and Zihao Wu (Kira 3); current developers Fabian Lange, Johann Usovitsch and
  Zihao Wu; see `AUTHORS`. References: P. Maierhöfer, J. Usovitsch, P. Uwer,
  Comput. Phys. Commun. 230 (2018) 99 [arXiv:1705.05610]; P. Maierhöfer,
  J. Usovitsch, Kira 1.2 release notes [arXiv:1812.01491]; J. Klappert,
  F. Lange, P. Maierhöfer, J. Usovitsch, Comput. Phys. Commun. 266 (2021)
  108024 [arXiv:2008.06494]; F. Lange, J. Usovitsch, Z. Wu, Comput. Phys.
  Commun. 322 (2026) 109999 [arXiv:2505.20197].
- **Base version:** Kira 3.1, upstream commit `5760f42` ("Set version number to 3.1").
- **License:** GPL-3.0-or-later; see the License section at the end of this file.
- **Delta versus upstream:** **+65/−19 lines of code across 5 files**
  (`src/kira/dataBase.cpp`, `src/kira/kira.cpp`, `src/kira/kira.h`,
  `src/kira/reduction.cpp`, `src/pyred/relations.h`), plus a one-line
  modification notice under the upstream copyright header of each of those
  files. The complete unified diff against upstream v3.1 is in
  [`PATCHES.diff`](PATCHES.diff).
- The FireFly subproject (J. Klappert, S. Y. Klein and F. Lange;
  arXiv:1904.00009, arXiv:2004.01463; <https://gitlab.com/firefly-library/firefly>,
  GPL-3.0-or-later) is **stock** (meson wrap `subprojects/firefly.wrap`,
  upstream branch `kira-2`; the builds validated here used FireFly commit
  `72e9e3c` of that branch). Nothing in FireFly is patched. `src/pyred`
  (touched by Patch 2) is Kira's own linear-algebra module pyRed; the Fermat
  pipe in `src/kira/connect2kira.*` derives from gateToFermat by Mikhail
  Tentyukov (GPL-2.0-or-later, relicensed from GPLv2-only with the author's
  agreement so that it could be included in Kira), as upstream's header states.

- **Not a patch, shipped alongside:** `bootloops-tools/` holds two GPL-3.0-or-later
  Python utilities written for the BootLoops toolkit by reading Kira's and
  FireFly's sources — `pyred_weight.py` (codec for pyRed's integral-weight
  layout, `integral_ordering` 1..8) and `ffsave_degree_census.py` (per-entry
  degree census of an `ff_save/` state, no solver work). They derive from
  `src/pyred/integrals.cpp`, `src/kira/ReadYamlFiles.cpp` and FireFly's
  save-state format, hence carry the GPL and the Kira/FireFly attribution in
  their headers; the Kira sources and build are unchanged by them. See
  `bootloops-tools/README.md`.

Neither patch has been submitted upstream yet. Both are small and we would be
glad to see them adopted; we intend to open a merge request with the Kira
developers once the deep-reduction speedup of Patch 2 is measured (see the
measurement note below), and this file will link it when it exists. We thank
the Kira developers for a code base in which a five-file, 65-line change was
all that our workload required.

## Patch 1 — SQLite transaction autocommit guards (`src/kira/dataBase.cpp`)

**What:** `DataBase::begin_transaction()` and `DataBase::commit_transaction()` now
consult `sqlite3_get_autocommit()` before issuing `BEGIN` / `COMMIT`.

- `begin_transaction()`: if a transaction is already open, the nested `BEGIN` is
  skipped (a second `BEGIN` on an already-open connection errors).
- `commit_transaction()`: `COMMIT` is only issued when a transaction is actually
  open. When SQLite reports autocommit there is nothing to commit, so the call
  becomes a no-op that just finalizes the prepared statement.

**Why:** several call sites issue begin/commit in tight loops, and an
intermediate `sqlite3_step()` error auto-rolls-back the transaction. During long
reductions this produced a "cannot commit - no transaction is active" error line
on essentially every commit — non-fatal, but repeated often enough in our
multi-day reductions to obscure other diagnostics. Both calls are now
idempotent.

## Patch 2 — parallel treatcoeff (`src/kira/kira.cpp`, `kira.h`, `reduction.cpp`, `src/pyred/relations.h`)

**What:** Kira's back-substitution coefficient normalization
(`Kira::treatcoeff2`, which pipes each coefficient through Fermat) is serialized
on a single mutex and Fermat instance (`fermat[0]`; simple and safe, and the
Fermat pool of size `coreNumber` that `initiate_fermat()` already spawns makes a
per-thread variant straightforward). The patch:

1. gives each pyred worker thread its own fixed Fermat pool slot
   (`thread_local` id assigned from an atomic counter, modulo the pool size);
2. rewrites the `RetrieveParallelMaster` treatcoeff cache wrapper in
   `relations.h` to probe the cache under the mutex, compute **outside** it, then
   insert-or-return-existing under the mutex; a cache hit holds the lock only for
   a map lookup rather than a Fermat round-trip. A cache miss may be computed by
   more than one thread, but Fermat canonicalization is deterministic, so
   duplicates are identical and the wasted work is bounded by the thread count
   per unique coefficient;
3. applies `--set_value` substitutions to **every** Fermat pool instance rather
   than only `fermat[0]`, since `treatcoeff2` now dispatches to `fermat[tid]`.

**Why:** on deep reductions the Fermat normalization phase serializes all worker
threads through one mutex and one external Fermat process; this removes that
bottleneck.

**Measurement note.** Validation: reduction outputs are
byte-identical to the stock binary across a 144-case comparison at 8 and 32
threads. The speedup, however, is **unproven at moderate depth**: on a benchmark
where the treatcoeff phase was only 7–17% of wall time (r = 12, s = 0 seeds),
the patched build ran at 0.98× stock, i.e. parity within noise. The patch
targets jobs deep enough (roughly r ≥ 13, s ≥ 1) for coefficient normalization
to dominate; measure before assuming a win.

## Build variants (no source delta — upstream meson options we rely on)

Standard build:

```sh
meson setup build --buildtype=release -Dfirefly=true -Dflint=true --prefix=<install>
ninja -C build install
```

- **jemalloc** (`-Djemalloc=true`): links the binary against jemalloc, which we
  found markedly more robust than glibc malloc for long, allocation-heavy
  reductions. `kira --version` reports "with jemalloc" when active. Requires the
  jemalloc development package. jemalloc is a documented Kira build option
  (`-Djemalloc=true`, default off); the Kira README reports "significantly
  increased performance, often by more than 20% from our experience if FireFly
  is used", with a caveat for some MPI setups.
- **kira128** (`-Dweight_width=128`): builds a 128-bit-integral-weight binary
  (meson names it `kira128`). This raises the weight ceiling from 64 to 128 bits,
  so the abort `Integral weight representation exceeds 64 bits`, which stock
  64-bit builds hit on very-high-`s` seed systems (observed at s = 437), no
  longer fires there. Validated by (a)
  running past the abort point on a job that kills the 64-bit binary, and (b) a
  regression job on which the 128-bit binary's reduction outputs (masters file
  and full reduction table) are byte-identical to the 64-bit binary's, at
  runtime parity (3.4 s vs 3.5 s). **Keep the 64-bit build as the default** —
  128-bit weights double the weight footprint for no benefit on normal jobs;
  reach for `kira128` only when you see the named abort.

Both options exist in upstream `meson_options.txt`; they are documented here
because this fork's production configuration is
`-Djemalloc=true` plus, for the affected job class only, `-Dweight_width=128`.

## Runtime dependency note

Kira's back-substitution requires **Fermat**, the computer algebra system by
Robert H. Lewis (Fordham University). Fermat is freeware distributed by its
author (the source of recent versions has been released under the GNU GPL on
request to him) and is **not** redistributed here — obtain it from
<http://home.bway.net/lewis/> and point `FERMATPATH` at the binary.

## License

Kira is free software, Copyright (C) 2017–2025 The Kira Developers (see
`AUTHORS`), distributed under the GNU General Public License, version 3 or (at
your option) any later version; the license text ships unmodified as `COPYING`.
This repository is a fork of Kira and is distributed under those same
GPL-3.0-or-later terms. The modifications described above are Copyright (c) 2026
Anthropic, PBC and are contributed under the GPL-3.0-or-later as well; they were
created by Matthew D. Schwartz, with the code written by Claude (Anthropic) under
his supervision. This is not an officially supported Anthropic product; it is
maintained by Matthew D. Schwartz (https://www.bootloops.ai). Each modified file carries a modification notice under its
upstream copyright header, as the GPL asks. Everything in this repository is under
Kira's GPL-3.0-or-later terms; the bundled gzstream (LGPL-2.1-or-later), SQLite
(public domain), the gateToFermat-derived Fermat pipe wrapper in
src/kira/connect2kira.* (M. Tentyukov, GPL-2.0-or-later) and the LaTeX/BibTeX
support files under doc/ (LPPL) keep their own terms exactly as upstream Kira
ships them. (The separate main BootLoops repository is MIT-licensed and contains
none of this code; it links here from upgrades/ENGINES.md.)
