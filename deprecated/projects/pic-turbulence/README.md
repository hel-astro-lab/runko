# PIC turbulence

This project uses the particle-in-cell (PIC) module to simulate driven kinetic plasma turbulence. The code uses solvers from `pic` and `emf` modules to time step the system.

A small test case of kinetic turbulence can be run with
```
mpirun -n 4 python3 pic.py --conf turb_small.ini
```

This uses the `pic.py` driver that loads a grid & problem configuration from `turb_small.ini`. 
The simulation uses 4 cores (`-n 4`).


You can visualize the box with
```
python3 plot3d.py --conf turb_small.ini --var rho --lap 200
python3 plot3d.py --conf turb_small.ini --var jz  --lap 200
```


The simulation can be analyzed by computing the energy evolution with
```
python3 plot_ene.py --conf turb_small.ini
```
that creates a plot into `output_dir/ene_histor.pdf`.


## Files

- `pic.py`: main driver script
- `init_problem.py`: parses configuration file and defines problem specific parameters
- `turb_small.ini`: simulation configuration file; a small turbulence setup that can be run on a laptop


and analysis files

- `plot3d.py`: plot emf fields on a 3d box
- `plot_ene.py`: computes and plots total energy in the system
- `units.py`: defines some standard units for the simulation analysis scripts

