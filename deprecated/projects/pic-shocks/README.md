# Collisionless shocks 

This project uses the particle-in-cell (PIC) module to simulate collisionless shocks. The code uses solvers from `pic` and `emf` modules to time step the system.

A test simulation can be run on laptop with
```
mpirun -n 4 python3 pic.py --conf 2dsig3.ini
```

This uses the `pic.py` driver that loads a grid & problem configuration from `2dsig3.ini`. 
The simulation uses 4 cores (`-n 4`).

The result can be analyzed with 
```
python3 plot_win_2d_shock.py --conf 2dsig3.ini --lap 1000
```
that visualizes the output file at a lap 1000 and saves it into the run folder.

To plot every frame, use 
```
./scripts.sh
```
which calls all the analysis scripts and loops over time slices.
Note that the script assumes zsh shell.


## Runko shock files

Simulation files consist of

- `pic.py`: main driver script
- `init_problem.py`: parses configuration file and defines problem specific parameters
- `2dsig3.ini`: simulation configuration file; a tiny shock setup that can be run on a laptop
- `injector.py`: python tool that controls the right wall and injects plasma on the fly
- `box_chunking.py`: python tool that controls the left wall and chunks away shock downstream to save computational resources
- `shock_toolset.py`: python io tool that tracks the shock front and saves a high-resolution field snapshots around the shock front

And analysis files of

- `plot_win_2d_shock.py`: analysis script that visualizes the shock.
- `plot_dens.py`: analysis script that plots the density as a ``mountain plot''
- `plot_jump_conds.py`: analysis script that calculates the shock jump conditions and plots them
- `plot_upstream_ene.py`: analysis script that calculates electromagnetic fields ahead of the shock
- `labellines.py`: utility function for plotting



