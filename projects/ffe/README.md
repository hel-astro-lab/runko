# FFE

Force-free electrodynamics code. Uses solvers from `ffe` and `em-fields` modules to time step the system.

A small test case can be run with

```
mpirun -n 4 python3 ffe.py --conf pkg_small.ini
```

This uses the `ffe.py` driver that loads a grid & problem configuration from `pkg_small.ini`. The simulation uses 4 cores (`-n 4`).

The simulation can be analyzed by computing the energy evolution with
```
python3 plot_ene.py --conf pkg_small.ini
```
and by rendering 3D snapshot (via `mayavi`) with
```
python3 plot_3d.py --conf pkg_small.ini --var bvec --lap 100
```
where
- `--var` defines the variable; `bvec` for magnetic field vector field
- `--lap` defines the simulation lap to be processed



## Files

- `ffe.py`: main driver script
- `problem.py`: parses configuration file and defines problem specific parameters
- `packet_setup.py`: defines `insert_em_fields` function that sets the initial electromagnetic field setup. Defines two colliding torsional Alfven wave packets.
- `units.py`: defines some standard units for the simulation analysis scripts
- `pkg_small.ini`: simulation configuration file; a small Alfven wave packet collision setup that can be run on a laptop
- `plot_ene.py`: computes and plots total energy in the system
- `plot_3d.py` renders 3D snapshots of the grid via `mayavi`


