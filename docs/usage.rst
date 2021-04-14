Code usage
############


Configuration files
===================

Simulation parameters are controlled by JSON configuration files.

Basic usage
-----------

The files have an `.ini` suffix and a simple structure of `parameter: value`. When fed into the drivers, Runko will automatically parse the file and construct a `Configuration` object, typically called `conf`. Any parameter in the configuration file gets turned into `conf.parameter` with a value of `value`; this makes it easy to add any new user-defined parameters to the driver scripts.

During the construction of `conf` the configuratoin file is also typically passed through a problem-spesific `init_problem.py` file that introduces a specialiced `Configuration_Problem` derived from the base class. It can do additional manipulation on the configuration file parameters before resulting in a `conf` object.


Typical parameters 
-----------

Typical default parameters found in configuration files include:

- [io]: basics
   - `outdir`: simulation output directory name (e.g., `run_1`)
   - `laprestart`: Simulation lap to restart from (`-1` = no restart, `0` = automatic, `X` lap to restart; `0` is default)
   - `restart`: No. of steps between restart files written; these restart files are cycled and overwritten during the simulation to save disk space.

- [io]: analysis
   - `interval`: No. of steps between analysis output files written
   - `stride`: thinning factor for analysis grids; results in saving only every `st`:th grid point in each dimension (typically either `1` for complete data snapshots or tile length for maximum compression)
   - `stride`: thinning factor for analysis grids; results in saving only every `st`:th grid point in each dimension (typically either `1` for complete data snapshots or tile length for maximum compression)

- [io]: deep io
   - `full_interval`: No. of steps between writing of full output files (these act as full restart files but are not overwritten; typically `-1` to save disk space) 

- [io]: shock io (optional and relevant only for shock simulations)
   - `box_nx`: length of a special grid that tracks the shock front during a simulation
   - `box_stride`: similar thinning factor as `stride` but for the special shock grid
   - `box_shift`: displacement of the box from the (automatically-located) shock front (use a value of `-box_nx/2` to center the grid on the shock front; `0` to track upstream).

- [simulation]: MPI control parameters (optional)
   - `mpi_track`: instead of partitioning the domain equally between processors we can optionally make a "caterpillar`-like track (along x-dimension) where tiles owned by different MPI ranks repeat themselves cyclically. This parameter controls the length of the track in units of tiles. A value of `128` would partition the tiles inside a sub-domain of `128 x Ny x Nz` equally between all the processors. The pattern is repeated for the length of the domain in x-direction. 
        - make sure there are reasonable number of tiles in one cycle per processor. Each rank has `mpi_track x Ny x Nz/Nprocs` continuous tiles in the subdomain, multiplied by the number of cycles, `Nx/mpi_track`.
          - The parameter can be tuned to roughly equal the active zone in a shock simulation. This way the simulation is always distributed among approximiately the same number of ranks. Physical length of the track is `mpi_track * NxMesh/c_omp`.



