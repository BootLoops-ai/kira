Kira 3.1
=============

New features
------------

 * Added the definition file `kira.def` to build an Apptainer container, targeted mainly at high-performance computing systems. See [merge request 51](https://gitlab.com/kira-pyred/kira/-/merge_requests/51).

Changes
-------

 * Changes to `numerical_points` (see [issue 91](https://gitlab.com/kira-pyred/kira/-/issues/91) and [merge request 50](https://gitlab.com/kira-pyred/kira/-/merge_requests/50)):
   - The results are now written out as Mathematica replacement rules.
   - Vanishing integrals are now written to the database and can be extracted with `kira2math`, `kira2form`, and `kira2formfill`.
   - The number of points for the output files can be specified with `numerical_points: {file: <file>, points_per_output_file: <points>}`. The old syntax `numerical_points: <file>` is still accepted.

 * Cleaned up and fixed the examples further.

Bug fixes
---------

 * Fixed that Kira crashed when an integralfamily was defined, but not reduced. See [issue 91](https://gitlab.com/kira-pyred/kira/-/issues/91).

 * Fixed that Kira could crash randomly when the option `zero_sectors` was set in `integralfamilies.yaml`. See [issue 93](https://gitlab.com/kira-pyred/kira/-/issues/93) and [merge request 48](https://gitlab.com/kira-pyred/kira/-/merge_requests/48).

 * Fixed that Kira crashed with the option `check_masters: true` if trivial integrals are in the list of preferred masters. See [issue 92](https://gitlab.com/kira-pyred/kira/-/issues/92) and [merge request 47](https://gitlab.com/kira-pyred/kira/-/merge_requests/47).

 * Fixed that Kira could not properly read CRLF files on Linux, resulting in unintended behavior or crashes. See [issue 88](https://gitlab.com/kira-pyred/kira/-/issues/88) and [merge request 53](https://gitlab.com/kira-pyred/kira/-/merge_requests/53).

 * Fixed a rogue comma in the output files produced with `numerical_points`. See [issue 91](https://gitlab.com/kira-pyred/kira/-/issues/91) and [merge request 50](https://gitlab.com/kira-pyred/kira/-/merge_requests/50).

 * Fixed that Kira crashed when the number of propagators provided in `integralfamilies.yaml` did not match its expectation and throw an error instead. See [merge request 52](https://gitlab.com/kira-pyred/kira/-/merge_requests/52).


Kira 3.0
=============

New features
------------

 * Added the job file option `truncate_sp:\n  - {topologies: [<topologies>], l: <n>}` where `<n>` can be chosen by the user. This improves Kira's seeding mechanics and is described in the Kira 3 publication in great detail. One can also find some guidelines on how to choose `<n>` in there.

 * Significantly improved the algorithm to select relevant equations for the reduction of the target integrals. The details can be found in the Kira 3 publication.

 * Allow users to add additional relations during the generation step with the job file option `extra_relations: <file>`, where `<file>` is the name of the file containing the additional relations in the format of user-defined systems. The additional relations are inserted after the relations for a specific sector were generated, or at the very end if they belong to sectors beyond the seeded sectors. An example is available in `examples/extra_relation`.

 * Allow the users to sample the system over phase space points in a prime field. When instructed with the job file option `numerical_points: <file>`, Kira will sample the system over the provided points and primes. The results are written to `results/<topology>/<file>_<prime>_<j>.m`. Kira writes at most 1000 entries per file and enumerates the files with the incremental number `<j>`. This feature is described in the Kira 3 publication in more detail and an example is available in `examples/numerical_IBP`.

 * The command line option `set_value` now affects the phase space point on which the system is generated. Furthermore, the user can also choose the prime with the command line option `prime_user=<prime>`.

 * The new job file option `symbols_2_nums: true` allows the user to write the generated system to disk with all variables replaced by the chosen numerical values set by `numerical_points` or `set_value`.

 * Kira can now be instructed to perform symbolic reductions. The option `symbolic_ibp: [<n1>, <n2>, ...]` in `config/integralsfamilies.yaml` can be used to treat the powers of the propagators `<ni>` as algebraic symbols. This feature is described in the Kira 3 publication in more detail and an example is available in `examples/symbolic_IBP`.

 * Added the job file option `check_masters: true/false` to check if the master integrals are a subset of those provided by the user with the options `select_masters` or `preferred_masters`. If not, Kira exits with an error message. See [issue 51](https://gitlab.com/kira-pyred/kira/-/issues/51) and [merge request 44](https://gitlab.com/kira-pyred/kira/-/merge_requests/44).

Changes
-------

 * Kira now requires a compiler with C++17 support.

 * Kira now generates IBP equations for sectors mapped away by symmetries. Compared to IBP equations, symmetry equations tend to grow much stronger in complexity with increasing values of `s`.

 * Kira now respects the seed specifications of `r`, `s`, and `d` set by the user for sectors related by symmetries and sectors containing preferred master integrals.

 * The examples have been restructured.

 * Kira now support 128-bit weights for its internal representation of integrals. This is required for projects with many dots, scalar products, or lines. This option requires 128-bit integer support from the compiler and can be set at compile time with `meson setup -Dweight_width=<width> build` where `<width>` can be set to 64 (default) or 128. For 128-bit, the executable is called `kira128`.

 * `permutation_option: 1` is now applied per default, but can manually be overridden by setting a different ordering scheme, 1-4, or setting `permutation: [<propagatornumber1>,<propagatornumber2>,...]`, see the changelog for version 2.3.

 * With the commandline option `-l` Kira now writes to `kira_DATE_TIME.log` instead of `kira.log`. See [issue 64](https://gitlab.com/kira-pyred/kira/-/issues/64) and [merge request 41](https://gitlab.com/kira-pyred/kira/-/merge_requests/41).

 * Changed how Kira handles errors and changed error messages in the process. See [issue 48](https://gitlab.com/kira-pyred/kira/-/issues/48) and [merge request 24](https://gitlab.com/kira-pyred/kira/-/merge_requests/24).

Bug fixes
---------

 * Fixed that Kira crashed with an out-of-bound access during the search for trivial sectors since version 2.3 when compiled with recent compiler versions. With older compiler versions it ran through, but some intermediate results were wrong, most notably the order of sectors written out to `sectormappings/nonTrivialSectors`. Final reduction tables should not have been affected. See [issue 78](https://gitlab.com/kira-pyred/kira/-/issues/78) and [merge request 42](https://gitlab.com/kira-pyred/kira/-/merge_requests/42).

 * Kira no longer freezes when communicating with Fermat 7.X. See [issue 75](https://gitlab.com/kira-pyred/kira/-/issues/75) and [merge request 46](https://gitlab.com/kira-pyred/kira/-/merge_requests/46). Thanks to Takahiro Ueda for the fix.

 * Fixed that vanishing coefficients in user-defined systems were not detected. See [issue 81](https://gitlab.com/kira-pyred/kira/-/issues/81).

 * Fixed that Kira did not complain if configuration files or files specified by the users were missing. See [issue 49](https://gitlab.com/kira-pyred/kira/-/issues/49) and [merge request 24](https://gitlab.com/kira-pyred/kira/-/merge_requests/24).

 * Fixed a rare deadlock which would lead to Kira being stuck forever. See [issue 69](https://gitlab.com/kira-pyred/kira/-/issues/69) and [merge request 38](https://gitlab.com/kira-pyred/kira/-/merge_requests/38).

 * Fixed that database connections were not properly closed. If many connections were established during a run, for example by reducing O(1000) topologies, the operating system then forbade Kira to open more files or pipes resulting in a crash. See [issue 74](https://gitlab.com/kira-pyred/kira/-/issues/74) and [merge request 39](https://gitlab.com/kira-pyred/kira/-/merge_requests/39).

 * Fixed a memory leak in Pak's algorithm when searching for symmetries. See [issue 79](https://gitlab.com/kira-pyred/kira/-/issues/79) and [merge request 43](https://gitlab.com/kira-pyred/kira/-/merge_requests/43).


Kira 2.3
=============

New features
------------

 * Added the option to reorder the propagators *internally* without requiring the user to redefine the integralfamily. The integrals in the output files will always stay in the order of the definition of the integralfamily. For some integralfamilies the difference in performance between good and bad orderings can be significant.
 The user can either choose a permutation manually or employ one of the four predefined ordering criteria. The manual order can be enforced by adding `permutation: [<propagatornumber1>,<propagatornumber2>,...]` to the family in `config/integralfamilies.yaml` where each `<propagatornumber>` is the number of a propagator from 0 to the number of propagators minus 1 *OR* 1 to the number of propagators, e.g. `permutation: [2,0,1]` corresponds to the internal order _propagator 2, propagator 0, propagator 1_. The automatic ordering can be enabled by adding `permutation_option: X` instead where `X=1,...,4`. In options 1 and 2 massless propagators are considered simpler than massive ones, in options 3 and 4 it is reversed. In options 1 and 3 propagators with less momenta will be considered simpler, in option 2 and 4 it is reversed.
 See [issue 34](https://gitlab.com/kira-pyred/kira/-/issues/34) and [merge request 28](https://gitlab.com/kira-pyred/kira/-/merge_requests/28).

 * Added the job option `run_initiate: input` to write the generated and selected system to `.kira` files in the directory `input_kira` to be used as a user-defined system. See [issue 43](https://gitlab.com/kira-pyred/kira/-/issues/43).

 * Added the option to use 128-bit weights so that Kira can reduce integrals with many dots and/or many scalar products. The number of sectors is still limited to 32 bits. Since there is a performance and memory impact, it has to be enabled with the compiler option `weight_width=128`. See [merge request 31](https://gitlab.com/kira-pyred/kira/-/merge_requests/31).

 * Extended the symmetry finder to work in some unusual cases.

Changes
-------

 * Improved the performance when exporting results. See [issue 31](https://gitlab.com/kira-pyred/kira/-/issues/31), [merge request 10](https://gitlab.com/kira-pyred/kira/-/merge_requests/10), and [merge request 27](https://gitlab.com/kira-pyred/kira/-/merge_requests/27).

 * The Autotools build system is deprecated, will no longer be supported, and will be removed in a future Kira release.

Bug fixes
---------

 * Fixed that user-defined systems with `otf: false` were not sorted after the selection run. For some highly unordered systems this could increase the runtime with `run_firefly` by orders of magnitude. See [issue 39](https://gitlab.com/kira-pyred/kira/-/issues/39) and [merge request 14](https://gitlab.com/kira-pyred/kira/-/merge_requests/14).

 * Fixed a segmentation fault for systems without `config` files. See [issue 41](https://gitlab.com/kira-pyred/kira/-/issues/41) and [merge request 13](https://gitlab.com/kira-pyred/kira/-/merge_requests/13).

 * Fixed that `insert_prefactors` did not work for user-defined systems with user-defined weights. See [issue 42](https://gitlab.com/kira-pyred/kira/-/issues/42) and [merge request 20](https://gitlab.com/kira-pyred/kira/-/merge_requests/20).

 * Fixed a segmentation fault when reducing a user-defined system with `run_firefly` and `iterative_reduction` if some requested integrals do not appear in the system at all. See [issue 44](https://gitlab.com/kira-pyred/kira/-/issues/44) and [merge request 16](https://gitlab.com/kira-pyred/kira/-/merge_requests/16).

 * Fixed that Kira did not export results if all selected integrals belong to trivial sectors. See [issue 66](https://gitlab.com/kira-pyred/kira/-/issues/66) and [merge request 34](https://gitlab.com/kira-pyred/kira/-/merge_requests/34).

 * Fixed that Kira did not export results for user-defined systems if all integrals vanish. See [issue 24](https://gitlab.com/kira-pyred/kira/-/issues/24) and [merge request 32](https://gitlab.com/kira-pyred/kira/-/merge_requests/32).

 * Killed some Fermat zombies. See [issue 45](https://gitlab.com/kira-pyred/kira/-/issues/45) and [merge request 19](https://gitlab.com/kira-pyred/kira/-/merge_requests/19).

 * Fixed that Kira exited with code 0 even if there was an error. See [issue 47](https://gitlab.com/kira-pyred/kira/-/issues/47) and [merge request 18](https://gitlab.com/kira-pyred/kira/-/merge_requests/18).

 * Fixed that `masters.final` was not written with `run_firefly`. See [issue 52](https://gitlab.com/kira-pyred/kira/-/issues/52) and [merge request 25](https://gitlab.com/kira-pyred/kira/-/merge_requests/25).


Kira 2.2
=============

New features
------------

 * Added the option `run_initiate: masters` to stop the reduction after the master integrals have been identified. See [issue 21](https://gitlab.com/kira-pyred/kira/-/issues/21).

 * The weight bits for user-defined systems with indexed integrals (T[a,b,c,...]) are now automatically adjusted. Before, Kira crashed if the weights did not fit into the default representation. See [issue 23](https://gitlab.com/kira-pyred/kira/-/issues/23).

Changes
-------

 * Changed the output format for the master integrals. This mainly affects the output to the console, `kira.log`, and `results/<TOPOLOGY>/masters(.final)`.

 * The command line option `--parallel`/`-p` now accepts the arguments `physical` and `logical` to exploit all physical or logical cores, respectively. It is still possible to manually choose any number of threads. If the option is given to Kira without argument, is uses `physical` by default. Under macOS, `physical` also uses the number of logical cores.

 * Added the missing file `master_equations` to `examples/master_equations/`. See [issue 33](https://gitlab.com/kira-pyred/kira/-/issues/33).

Bug fixes
---------

 * Fixed that Kira 2.1 missed some symmetries and found too many master integrals. See [issue 22](https://gitlab.com/kira-pyred/kira/-/issues/22).

 * Fixed that `iterative_reduction` could select wrong equations resulting in long runtimes and/or only partially reduced results. See [issue 25](https://gitlab.com/kira-pyred/kira/-/issues/25).

 * Fixed that the result of `run_back_substitution` was exported wrongly with `kira2form` and `kira2formfill`. See [issue 28](https://gitlab.com/kira-pyred/kira/-/issues/28).

 * Check whether the file for the options `select_masters` and `preferred_masters` can actually be found instead of silently ignoring the file. See [issue 33](https://gitlab.com/kira-pyred/kira/-/issues/33).

 * Fixed a crash with `Fermat error: \n*** Inappropriate symbol` in `run_triangular` and `run_back_substitution`. See [issue 29](https://gitlab.com/kira-pyred/kira/-/issues/29).

 * Fixed that Kira might stop a reduction with multiple sectors with the error message `Seed specification: rmax must not be lower than the number of lines of the sector.`. See [issue 17](https://gitlab.com/kira-pyred/kira/-/issues/17).

 * Fixed that Kira did not terminate in a reduction of an user-defined system if there are no master integrals. See [issue 24](https://gitlab.com/kira-pyred/kira/-/issues/24).

 * Fixed that it was not always possible to resume a reduction of a user-defined system with multiple topologies. See discussion in [issue 30](https://gitlab.com/kira-pyred/kira/-/issues/30).


Kira 2.1
=============

New features
------------

 * Added `kira2formfill` to export `fill` statements to `results/<TOPOLOGY>/kira_<SELECTED>.h`, which can be used to fill `TableBase` objects in FORM.

Changes
-------

 * Introduced a version number for the database. Databases produced prior to Kira 2.1 cannot be used for `kira2form` and `kira2formfill` unless they are upgraded with the command line argument `--force_database_format=fermat|firefly`. Read this [wiki article](https://gitlab.com/kira-pyred/kira/-/wikis/Upgrading-the-database-format) for details.

 * The numerators of the rational functions exported with `kira2form` and `kira2formfill` are now enclosed by `num()`.

 * The prefactors submitted through `insert_prefactors` are now enclosed by `prefactor[]` when employing `kira2math` and by `prefactor()` when using `kira2form` and `kira2formfill`.

 * `Meson` build system: Removed the default option `libdir=lib`. On Debian-based systems (e.g. Ubuntu) this should be set manually by the user as explained in [`README.rst`](https://gitlab.com/kira-pyred/kira/-/blob/master/README.rst).

Bug fixes
---------

 * In a reduction with Fermat, integers and rational numbers were truncated after 2048 characters leading to wrong results or a segmentation fault if the 2048th character was the `/` character. This was fixed. See [issue 14](https://gitlab.com/kira-pyred/kira/-/issues/14).

 * Fixed that using different loop momenta in different topologies led to wrong symmetries. See [issue 7](https://gitlab.com/kira-pyred/kira/-/issues/7).

 * Fixed that loop momenta defined by appending a `2` to another momentum, e.g. `k` and `k2`, led to wrong trivial sectors and symmetries. See [issue 8](https://gitlab.com/kira-pyred/kira/-/issues/8).

 * Fixed that in very rare cases wrong symmetries were found when defining propagators in the Euclidean metric. See [issue 10](https://gitlab.com/kira-pyred/kira/-/issues/10).

 * Fixed that exporting results produced with `run_firefly` employing `kira2form` could produce wrong `id` statements. See [issue 11](https://gitlab.com/kira-pyred/kira/-/issues/11).

 * Fixed that some files generated by `run_triangular` were not properly read in by `run_firefly: back` which then only solved a subsystem and produced an incomplete result. See [issue 9](https://gitlab.com/kira-pyred/kira/-/issues/9).

 * Fixed that `kira2file` added extra lines consisting of `0` after vanishing integrals so that the resulting files could not be read in by `input_system`. See [issue 15](https://gitlab.com/kira-pyred/kira/-/issues/15).

 * Fixed that pure export jobs could crash with `pyred::init_error`. See [issue 16](https://gitlab.com/kira-pyred/kira/-/issues/16).

 * Fixed that Kira generates and reduces the **whole** system of equations for **all** integrals within the seed boundaries if all selected integrals belong to trivial sectors or are master integrals. See [issue 12](https://gitlab.com/kira-pyred/kira/-/issues/12).

 * `results/<TOPOLOGY>/masters` is now written in its entirety before generating the `tmp/<TOPOLOGY>/SYS` files.
