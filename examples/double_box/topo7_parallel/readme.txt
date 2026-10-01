
Example how to parallelise reduction across several machines
============================================================

Run the forward elimination. Note that --parallel=n enables parallelisation on
one machine (number of CPU cores of the machine).

$ cd examples/topo7_parallel
$ kira --parallel=4 prepare.yaml

When the forward elimination finished, look into the file results/topo/masters.
We already created job files to run on two machines in parallel. In each job
file (parallel1.yaml and parallel2.yaml), enter the numbers (starting with 1) of
the master integrals for which the reduction should be done under the option
"select_masters_reduction" (for demonstration purposes this has already been
done in this example; parallel1.yaml: odd numbers, parallel2.yaml: even numbers
up to 31 which is the total number of master integrals in this example).

Copy the entire folder topo7_parallel to examples/topo7_parallel2.

$ cd ..
$ cp -r topo7_parallel topo7_parallel2

Run parallel1.yaml (on machine 1) in topo7_parallel.

$ cd examples/topo7_parallel
$ kira --parallel=4 parallel1.yaml

Run parallel2.yaml (on machine 2) in topo7_parallel2.

$ cd examples/topo7_parallel2
$ kira --parallel=4 parallel2.yaml

When both jobs finished, merge the databases.

$ cd examples/topo7_parallel
$ kira merge.yaml

This merges the database topo7_parallel2/results/kira.db into
topo7_parallel/results/kira.db.

You will get the following warning messages from SQLite which can be ignored:
  SQLite3: table EQUATION already exists
  SQLite3: table INTEGRALORDERING already exists

Afterwards, topo7_parallel2 can be deleted.

From the merged database, integrals can be extracted with kira2math or kira2form
as usual, e.g. with

$ kira export.yaml
