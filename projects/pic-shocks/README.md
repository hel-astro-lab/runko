# Collisionless shocks 

This project uses the particle-in-cell (PIC) module to simulate collisionless shocks. The code uses solvers from `pic` and `em-fields` modules to time step the system.

A small test simulation can be run with
```
mpirun -n 4 python3 pic.py --conf shock_mini.ini
```

This uses the `pic.py` driver that loads a grid & problem configuration from `shock_mini.ini`. 
The simulation uses 4 cores (`-n 4`).

The simulation result can be analyzed by 
```
python3 prtcl_spec.py --conf shock_mini.ini
```
that computes the particle spectra and saves the output into `out-mini/prtcl_spec.pdf`.

To plot every frame, use
```
python3 plot_shock.py --conf shock_mini.ini
```


## Files

- `pic.py`: main driver script
- `problem.py`: parses configuration file and defines problem specific parameters
- `shock_mini.ini`: simulation configuration file; a tiny shock setup that can be run on a laptop
- `shock_small.ini`: simulation configuration file; a small shock setup starts to show real plasma dynamics
- `prtcl_spec.py`: analysis script that computes particle spectra.
- `plot_shock.py`: analysis script that visualizes the 2D fields.


