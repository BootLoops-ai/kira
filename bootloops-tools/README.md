# bootloops-tools — two GPL Python utilities that read Kira's and FireFly's files

This directory holds two small Python programs written for the BootLoops
toolkit that decode files produced by Kira and FireFly. They live in this
repository, and not in the MIT-licensed `bootloops` toolkit, because they
were written by reading Kira's and FireFly's GPL sources: the integral-weight
layout follows `src/pyred/integrals.cpp` (`Topology::sector_weight_table`,
`Integral::to_weight`, `Integral::weight_to_dots`) and the `integral_ordering`
mapping of `src/kira/ReadYamlFiles.cpp`, and the `ff_save` reader follows the
on-disk state format that FireFly's `RatReconst::save_state` writes. They are
therefore derived works and are distributed under the **GNU General Public
License, version 3 or later** (`../COPYING`), the license of Kira and FireFly,
with attribution to the Kira Developers (Copyright (C) 2017-2025, see
`../AUTHORS`) and to the FireFly authors (J. Klappert, S. Y. Klein, F. Lange).
The Python port and its extensions are Copyright (c) 2026 Anthropic, PBC; they
were created by Matthew D. Schwartz, with the code written by Claude (Anthropic)
under his supervision. Neither file is part of Kira itself and the Kira build
does not use them.

| file | what it does |
|---|---|
| `pyred_weight.py` | Codec between Kira integral names (`fam[a1,...,aN]` index tuples) and the integer weights Kira stores in `results/kira.db` and hands to FireFly, for the full `integral_ordering` range 1..8: sector rank (both sector orderings, out-of-tree sectors ranked after in-tree ones), the four dots/scalar-product field layouts, and the composition indices, with the bit widths read from the `WEIGHTBITS` row of `kira.db`. Encoding asserts every field fits; decoding fails loud (`ValueError`) on any inconsistent weight. A library module (`WeightCodec`, `from_kiradb`), no command line. |
| `ffsave_degree_census.py` | Reads the `ff_save/` state of a running, stopped or never-finishing Kira + FireFly reduction and writes a per-entry census, (target, master) -> (max degree of numerator, max degree of denominator, done flag), with per-variable degrees for multivariate runs, without any new solver work. Entries are named through the weight codec (file maps first, inverse decode otherwise); a round-trip layout gate (`decode(encode(name)) == name` on every known name) runs before anything is reported. |

Usage:

```sh
# degree census of a reduction directory (snapshot ff_save/states first if Kira is still running)
python3 bootloops-tools/ffsave_degree_census.py <rundir> --family <fam> [--vars d,s,t] [--workers 16] [--out census.json]

# the codec as a library
python3 -c "import sys; sys.path.insert(0, 'bootloops-tools'); import pyred_weight as pw; \
            c = pw.WeightCodec(4, (8, 8, 12, 12), 8); print(c.encode((2, 1, 0, -1)))"
```

`<rundir>` must contain `ff_save/states/`, the target list (`de_targets` or
`target`), `preferred`, `results/<family>/masters` and `results/kira.db`.
Both scripts need only the Python 3 standard library.

**How the MIT toolkit finds them.** The `kira-stack` and `ffcapital` packages
of `bootloops` use these two files optionally, through the small loader
`tools/kira-stack/kira_gpl_tools.py` there. The loader looks in the directory
named by the environment variable `BOOTLOOPS_KIRA_TOOLS` if it is set, and
otherwise in `../kira/bootloops-tools` relative to the
`bootloops` checkout, which is where this directory sits when the two
repositories are checked out side by side. When neither is present the
toolkit's battery legs that need them skip by name and the routes that need
them stop with an `ImportError` naming this repository; nothing falls back
silently.

```sh
# side-by-side checkouts: nothing to set
ls ../kira/bootloops-tools/pyred_weight.py
# any other layout
export BOOTLOOPS_KIRA_TOOLS=/path/to/kira/bootloops-tools
```

## License

GPL-3.0-or-later, as stated in each file's header; see `../COPYING`. Please
credit and cite Kira and FireFly (references in `../README.rst`) when these
utilities contribute to published work.
