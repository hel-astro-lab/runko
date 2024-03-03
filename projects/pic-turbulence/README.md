# PIC turbulence

This project uses the particle-in-cell (PIC) module to simulate kinetic turbulence. The code uses solvers from `pic` and `em-fields` modules to time step the system.

A small test case of kinetic turbulence can be run with
```
mpirun -n 4 python3 pic.py --conf turb_small.ini
```

This uses the `pic.py` driver that loads a grid & problem configuration from `turb_small.ini`. 
The simulation uses 4 cores (`-n 4`).

The simulation can be analyzed by computing the energy evolution with
```
python3 plot_ene.py --conf turb_small.ini
```
that creates a plot into `turb1/ene_histor.pdf`.


## Files

- `pic.py`: main driver script
- `problem.py`: parses configuration file and defines problem specific parameters
- `units.py`: defines some standard units for the simulation analysis scripts
- `turb_small.ini`: simulation configuration file; a small turbulence setup that can be run on a laptop
- `plot_ene.py`: computes and plots total energy in the system

