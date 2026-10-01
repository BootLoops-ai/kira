
==========================
 Kira 1.1 - Release Notes
==========================

---------------
 26 April 2018
---------------

Bug fixes
=========

Occasional segmentation faults
------------------------------

In specific circumstances (e.g. running the example ``tennisCourt``) Kira
crashed because of a segmentation fault. This was due to an internal limit
affecting the memory management which was hard coded using a preprocessor macro.
In the new Kira version this limit is now adjusted automatically during runtime.
We thank Pascal Wasser for pointing this problem out to us.

Note that no wrong results were produced. Results obtained with the previous
version are thus not affected.

Minor changes
=============

Number of Fermat workers
------------------------

In the previous Kira version the maximal number of Fermat workers was 63. In the
new version the number of Fermat workers can now be freely chosen by the user
and is practically only limited by the available main memory.

Path to the Fermat executable
-----------------------------

The program now gives an error message if the environment variable
``FERMATPATH`` points to a directory rather than to the ``Fermat`` executable
itself. In the previous version this lead to a segmentation fault.

kira.db moved
-------------

To use symmetry relations among different topologies, the previous Kira version created a copy of the ``kira.db`` file for each connected topology. To save disk space the new version now writes all results to one common file ``kira.db`` located in the directory ``results``. If in addition the option ``data_file: false`` is set in the job file, Kira writes only the database file and skips writing the results in human readable form (e.g. files ``kira`` and ``id2int``).

Performance issue
-----------------

In the old version a performance problem was encountered when coefficients of
serval GB in size occured during the reduction. This was due to an inefficient
allocation of new memory and has been fixed in version 1.1.


New Features
============

Usage examples for all new features can be found in the example directory
``examples/topo7`` which is included in the package. The new features are
illustrated in ``jobs4.yaml``.

Advanced seed selection
-----------------------

A new option for the seed selection of propagator powers has been introduced.
The seed specification ``r: [t,t+n]`` selects all seed integrals with at most
``n`` dots (``r-t``) on each sector, independently of the number of lines in the
sector.

This can be used to generate a smaller system of equations than with the
``r: [t,rmax]`` option. Note that if the new option is used in conjunction with
``select_integrals``, it will usually not affect the reduction time
significantly unless the equation generator is a bottleneck.

Example::

   jobs:
     - reduce_sectors:
         sector_selection:
           select_recursively:
             - [topo7,127]
         identities:
           ibp:
             - {r: [t,7], s: [0,1]}
             - {r: [t,t+2], s: [0,1]} # New feature

Cut propagators
---------------

The new option::

   cut_propagators: [n1,n2,n3,...]

where ``[n1,n2,n3,...]`` are the numbers of the cut propagators has been
introduced in the topology definition file ``integralfamilies.yaml``. Propagator
numbers start with 1.

Example::

   integralfamilies:
     - name: "topo7"
       loop_momenta: [k1,k2]
       top_level_sectors: [127]
       propagators:
         - ["-k1", 0]         #1
         - ["k2", 0]          #2
         - ["-k1+k2", 0]      #3
         - ["k1+q2", "m2^2"]  #4
         - ["k2-p2", 0]       #5
         - ["-k1+p1+p2", 0]   #6
         - ["k2-p1-p2", 0]    #7
         - ["k1-p2", 0]       #8
         - ["k2-q2", 0]       #9
       cut_propagators: [3,4]

Here, the 3rd and 4th propagators will be treated as cut propagators. This means
that during the reduction all integrals in which a cut propagator has
non-positive power are set to zero.

Basis change
------------

By default, Kira chooses the set of master integrals automatically according to
its internal integral ordering. Since version 1.1 it is now possible to specify
a list of integrals which should preferably be used as master integrals. The
list should be written to a file in the following format::

  # file name: basisChange
  - [topo7,0,0,2,1,1,0,0,0,0]
  - [topo7,0,2,1,1,0,1,0,0,0]
  - [topo7,2,0,1,1,1,1,0,0,0]
  ...

The file should be located in the job directory. To make use of it, set the
following option in the job file::

   select_integrals:
     select_masters: "basisChange"

Note that if additional master integrals are required, they will be chosen
automatically. Furthermore, integrals from the list which are linearly dependent
will be eliminated. In the current version no mechanism is forseen to change the
basis of master integrals after the reduction.

Anchor points
-------------

The results of the backward substitution will be commited to the database
``kira.db`` every 10 minutes. With version 1.1 the backward substitution can be
stopped at any time. To resume the backward substitution at a later time and to
load the results from a previous run the following option must be set in the job
file::

   conditional: true

An option to interrupt the backward substitution gracefully will be provided in
a future Kira version.

Calculate reduction coefficients only for selected master integrals
-------------------------------------------------------------------

This feature allows Kira to focus the reduction on a subset of the master
integrals. The coefficients of this subset of master integrals will be
determined, whereas the coefficients of all other master integrals will be
treated as zero. To perform the complete reduction, it must be run several times
until the coefficients of all master integrals are calculated. To make use of
this feature, set the following option in the job file::

   select_integrals:
     select_masters_reduction:
       - [topo7, [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31]]
       - [topo7, [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]]

This option only affects the back substitution and it must be used in
combination with the option ``select_mandatory_list`` or
``select_mandatory_recursively``. When the reduction process ``run_pyred``
terminates, the list of master integrals will be written to the file ``masters``
in the ``results`` directory. During the ``back_substitution`` reduction step
this file will be read and all master integrals will be enumerated from 1 to N,
where N is the number of master integrals. With the option
``select_masters_reduction`` one may give a list of integers from 1 to N (a
subset of {1,...,N}) which refers to the master integrals in the file
``masters``. In the example above we give two lists of integers, which refer to
the complete set of 31 master integrals of the topology ``topo7``. Here, Kira
will perform the reduction to the master integrals with odd numbers first and
subsequently to master integrals with even numbers. There a two applications of
this option:

* Reduce the memory consumption during the back substitution by calculating
  coefficients for different master integrals sequentially.
* Parallelise the reduction across different machines by letting each machine
  handle the coefficients of a subset of the master integrals.

Note that, at present, ``select_masters_reduction`` cannot be used in
conjunction with the option ``conditional: true`` [#]_.

.. [#] The reason is that, on resume, Kira needs to determine which integrals are still unreduced, but currently it cannot distinguish between integrals which are reduced to a complete set of master integrals and integrals in which some of the coefficients are set to zero. A workaround to make the option ``conditional:true`` usable in this case is to rename the file ``kira.db`` before the next reduction to a new subset of master integrals. The user will then need to merge the reduction formulas from the different database files manually.
