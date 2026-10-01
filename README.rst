**BootLoops fork.** This is the BootLoops project's lightly patched fork of
`Kira <https://gitlab.com/kira-pyred/kira>`__, the Feynman-integral reduction
program created by Philipp Maierhöfer, Johann Usovitsch and Peter Uwer and
developed today by Fabian Lange, Johann Usovitsch and Zihao Wu, with Jonas
Klappert (Kira 2.0); Philipp Maierhöfer contributed through Kira 2.3 and Peter
Uwer through Kira 1.1 (see ``AUTHORS``). The fork is a snapshot at upstream
release tag ``3.1`` (commit ``5760f42f8decd7eb5474c780db75b5e5e7b450a0``,
"Set version number to 3.1") with our modifications applied; upstream's
history and releases live at https://gitlab.com/kira-pyred/kira. What changed
and why: `PATCHES.md <PATCHES.md>`__ (the complete unified diff against the
base is `PATCHES.diff <PATCHES.diff>`__; +65/-19 lines of code in five source
files). ``bootloops-tools/``: two GPL Python utilities derived from Kira's
integral-weight layout and FireFly's save format (a weight codec and an
``ff_save`` degree census, used optionally by the ``bootloops`` toolkit;
see `bootloops-tools/README.md <bootloops-tools/README.md>`__); they are not
part of Kira and the build does not touch them. Upstream's license
(GPL-3.0-or-later, ``COPYING``) and notices (``AUTHORS``, copyright headers)
are unchanged; our modifications are released
under the same GPL-3.0-or-later license, and each modified file carries a
modification notice under its upstream copyright header. The modifications are
Copyright (c) 2026 Anthropic, PBC; created by Matthew D. Schwartz, with the code
written by Claude (Anthropic) under his supervision. This is
not an officially supported Anthropic product; it is maintained by Matthew D. Schwartz
(https://www.bootloops.ai). (This fork is part of the BootLoops release; the main
BootLoops repository, ``bootloops``, published beside it by the same
organization, is MIT-licensed and does not include this fork; it points here from
``upgrades/ENGINES.md``. In this repository everything is under Kira's own
GPL-3.0-or-later terms, with the bundled gzstream and SQLite sources keeping
the terms upstream ships them under.)

This fork is not affiliated with or endorsed by the Kira developers, to whom
we are grateful: Kira with FireFly carried a large share of the reductions in
the BootLoops program, including every reduction inside our AMFlow.cpp fork.
Please credit and cite Kira itself when using this fork: P. Maierhöfer,
J. Usovitsch, P. Uwer, *Kira — a Feynman integral reduction program*,
Comput. Phys. Commun. 230 (2018) 99 [arXiv:1705.05610]; J. Klappert, F. Lange,
P. Maierhöfer, J. Usovitsch, *Integral reduction with Kira 2.0 and finite
field methods*, Comput. Phys. Commun. 266 (2021) 108024 [arXiv:2008.06494];
F. Lange, J. Usovitsch, Z. Wu, *Kira 3: integral reduction with efficient
seeding and optimized equation selection*, Comput. Phys. Commun. 322 (2026)
109999 [arXiv:2505.20197]; and, when the finite-field mode is used, FireFly:
J. Klappert, F. Lange, *Reconstructing rational functions with FireFly*,
Comput. Phys. Commun. 247 (2020) 106951 [arXiv:1904.00009]; J. Klappert,
S. Y. Klein, F. Lange, *Interpolation of dense and sparse rational functions
and other improvements in FireFly*, Comput. Phys. Commun. 264 (2021) 107968
[arXiv:2004.01463]. Kira's algebraic back-substitution uses Fermat, the
computer algebra system by Robert H. Lewis of Fordham University
(http://home.bway.net/lewis/), which is not distributed here. Report problems
with this fork on this repository's issue tracker, not to the Kira developers.
The upstream README follows unchanged.

-----------------------------------
Maintenance, reporting and security
-----------------------------------

**Maintenance.** This repository is maintained by Matthew D. Schwartz, not by
Anthropic. It is not an officially supported Anthropic product, and Anthropic
does not provide support, updates or fixes for it.

**Reporting issues.** Please report bugs and security problems through this
repository's GitHub issues.

**Security considerations.** Treat input files from others as code. These are
research tools meant to be run locally on inputs you trust. Many of them
evaluate the contents of their input files (JSON, YAML, ``.m``, ``.ms``, ``.jl``,
pickle and similar), so a file received from someone else can run arbitrary
commands on your machine. Only run files you wrote yourself or got from a
source you trust, or run them in a sandbox or container. The integrity checks
and certificates in this repository guard against accidents. They are not a
security boundary.


.. sectnum::

=============================================
 Kira - A Feynman Integral Reduction Program
=============================================

.. --------
..  Readme
.. --------

.. contents:: Table of Contents
..    :depth: 2

.. Release notes
.. =============
..
.. See `ChangeLog <https://gitlab.com/kira-pyred/kira/blob/master/ChangeLog>`_.

Installation
============

Statically linked executables of ``Kira`` releases for ``Linux x86_64`` are available on https://kira.hepforge.org. This executable has all optional features included except for ``MPI``. If you require ``MPI``, you must compile ``Kira`` yourself against the ``MPI`` version used on your computer cluster.

To build ``Kira`` from source, clone the ``release`` branch of the ``Git`` repository with

::

   git clone https://gitlab.com/kira-pyred/kira.git -b release

to obtain the latest release version, or the ``master`` branch with

::

   git clone https://gitlab.com/kira-pyred/kira.git

to obtain the lastest development snapshot. Specific releases are available as ``Git tags``, e.g. ``kira-2.0``.

Prerequisites
-------------

**Platform requirements**

``Linux x86_64`` or ``macOS``.

**Compiler requirements**

A ``C++`` compiler supporting the ``C++17`` standard and a ``C`` compiler supporting the ``C11`` standard.

**Build system requirements**

``Kira`` can either be built with the ``Meson`` build system (http://mesonbuild.com) version ``0.46`` (or later) and ``Ninja`` (https://ninja-build.org), or with the ``Autotools`` build system. We recommend to use the ``Meson`` build system.

If ``Meson`` is not available on your system, it can be installed with

::

   pip3 install --user meson

(requires ``Python 3.5`` or later). With the option ``--user``, ``Meson`` can be installed as unprivileged user into the home directory (``~/.local/bin``). The ``Ninja`` binary can be downloaded from https://ninja-build.org.

**Dependencies**

``Kira`` requires the following packages to be installed on the system:

- ``GiNaC`` (https://www.ginac.de), which itself requires ``CLN`` (https://www.ginac.de/CLN),
- ``zlib`` (http://zlib.net).
- ``Fermat`` (http://home.bway.net/lewis) is required to run ``Kira``.

If the ``Fermat`` executable is not found automatically at startup, or a specific ``Fermat`` installation should be used, the path to the ``Fermat`` executable can be provided via the environment variable ``FERMATPATH``.

Depending on the enabled optional features of ``Kira``, the following packages are required in addition:

- ``GMP`` (https://gmplib.org) if ``FireFly`` is used,
- ``MPFR`` (https://www.mpfr.org) if ``FLINT`` is used,
- ``MPI`` (disabled by default) for parallelization on computer clusters,
- ``jemalloc`` (disabled by default) (http://jemalloc.net) for more efficient memory allocation.

The following dependencies can be automatically built and installed as subprojects with the ``Meson`` build system. I.e. if they are not found on the system, they will be built automatically along with ``Kira``.

- ``yaml-cpp`` (required) (https://github.com/jbeder/yaml-cpp),
- ``FireFly`` (optional, enabled by default). If you decide to install ``FireFly`` manually, we recommend to use the version from the branch ``kira-2`` of the ``Git`` repository at https://gitlab.com/firefly-library/firefly. This branch will remain compatible with ``Kira 2``.
- ``FLINT`` (optional, enabled by default) (http://www.flintlib.org). We recommend using ``FLINT``, because it not only offers better performance for the finite field arithmetic, but is also required to enable some features of ``FireFly``, most notably the factor scan.

If the ``Autotools`` build system is used, all enabled dependencies must be installed manually. If ``FireFly`` is not built as a subproject, to use ``FLINT`` and ``MPI``, they must be enabled in ``FireFly``'s ``CMake`` build system.

Note that ``GiNaC``, ``CLN``, ``yaml-cpp``, and ``FireFly`` must have been compiled with the same compiler which is used to compile ``Kira``. Otherwise the linking step will most likely fail. If you are using the system compiler, you can usually install ``GiNaC``, ``CLN``, and ``yaml-cpp`` via your system's package manager. However, if you are using a different compiler, in practice this usually means that you also have to build these packages from source and, if installed with a non-default installation prefix, the environment variables ``CPATH``, ``LD_LIBRARY_PATH`` and ``PKG_CONFIG_PATH`` must be set accordingly.

Compiling Kira with the Meson build system
------------------------------------------

To build ``Kira`` with the ``Meson`` build system, ``Meson 0.46`` (or later) and ``Ninja`` are required.

To compile and install ``Kira``, run

::

   meson setup --prefix=/install/path builddir
   cd builddir
   ninja
   ninja install

where ``builddir`` is the build directory. Specifying the installation prefix with ``--prefix`` is optional. On Debian-based systems (e.g. Ubuntu) one may want to use ``--libdir=lib``, otherwise libraries will be installed in ``$PREFIX/lib/<archdir>`` where ``<archdir>`` is e.g. ``x86_64-linux-gnu``.

**Build options**

- ``-Dfirefly=false`` (default: ``true``): Build without ``FireFly`` support.
- ``-Dflint=false`` (default: ``true``): If ``FireFly`` is built as a subproject, disable ``FLINT``.
- ``-Dmpi=true`` (default: ``false``): If ``FireFly`` is built as a subproject, enable ``MPI``. This is known to work best with ``OpenMPI`` (https://www.open-mpi.org). For performance reasons, we recommend ``MPICH`` (https://www.mpich.org), though.
- ``-Dcustom-mpi=<name>``: If your ``MPI`` installation provides a ``pkg-config`` file, but is not found automatically with ``-Dmpi=true``, pass the name of the ``MPI`` implementation as ``<name>``, e.g. ``-Dcustom-mpi=mpich``. Some systems don't provide a ``pkg-config`` file for ``MPICH``. In that case we recommend to install ``FireFly`` with its own ``CMake`` build system instead.
- ``-Djemalloc=true`` (default: ``false``): Link with the ``jemalloc`` memory allocator. This can lead to significantly increased performance, often by more than 20% from our experience if ``FireFly`` is used. However, using ``jemalloc`` may not work on some systems, especially in combination with certain ``MPI`` implementations (This can depend on subtleties like the linking order of the ``jemalloc`` and ``MPI`` libraries). Alternatively, to use ``jemalloc``, one can set the environment variable ``LD_PRELOAD`` to point to ``jemalloc.so`` and export it.
- ``-Dweight_width=128`` (default: ``64``): Build with 128-bit integers for the internal weight representation of integrals. This is required for projects with many dots, scalar products, or lines. This option requires 128-bit integer support from the compiler. For 128-bit, the executable is called ``kira128``.

To show the full list of available build options, run
``meson configure`` in the build directory.

**Subprojects**

If ``yaml-cpp`` or ``FireFly`` are not found on the system, per default they will be downloaded and built as ``Meson`` subprojects. If the option ``-Dflint=true`` (default) is set and ``FireFly`` is built as a subproject, also ``FLINT`` will be downloaded and built as a subproject if it is not found on the system.

The usage of subprojects can be controlled with the following options:

- ``--wrap-mode=nodownload``: Do not download subprojects, but build them if already available (and not found on the system).
- ``--wrap-mode=nofallback``: Do not build subprojects, even if the libraries are not found on the system.
- ``--wrap-mode=forcefallback``: Build subprojects even if the libraries can be found on the system.
- ``--force-fallback-for=<deps>``: Like ``forcefallback``, but only for dependencies in the comma separated list ``<deps>``. Overrides ``nofallback`` and ``forcefallback``.

These options are only fully supported with ``Meson 0.49`` or later. For details see https://mesonbuild.com/Subprojects.html.

*Note:* Subprojects are not updated automatically. To update subprojects, run

::

   meson subprojects update

(requires ``Meson 0.49`` or later). ``Git`` subprojects can of course also be manually updated by running

::

   git pull

in the corresponding subproject directory (e.g. ``subprojects/firefly``).

Compiling Kira with the Autotools build system
----------------------------------------------

*Note:* The Autotools build system is deprecated and will be removed in a future Kira release. Please use the Meson build system.

To compile and install ``Kira`` with the ``Autotools`` build system, first run

::

   autoreconf -i

and then compile and install with

::

   ./configure --prefix=/install/path --enable-firefly=yes
   make
   make install

where the optional ``--prefix`` argument sets the installation prefix. Without the option ``--enable-firefly=yes``, ``Kira`` will be built without ``FireFly`` support. Note that subproject installation is not supported with the ``Autotools`` build system, i.e. all dependencies must be installed manually.
