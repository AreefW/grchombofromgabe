## What is this example?
This example is a chombo code to rebuild checkpoint file from `GABRel` which has one big box to `GRChombo` checkpoint file which has many boxes. Doing this will make the parallelisation in `GRChombo` being done correctly
## To run
```mpirun -np {8} ./rebuild_checkpoint old_checkpoint.h5 new_checkpoint.3d.hdf5 {min_box_size} {max_box_size}```