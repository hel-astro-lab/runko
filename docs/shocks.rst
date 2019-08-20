Running the shock and tools
############

In this section it is explained how to run Runko and the given analytical tools

Running a first shock
===================
You will find the scripts to run the simulation and the analytical tools in

.. code-block:: bash

   /runko/projects/shocks
   
You can use this directory or make your own using

.. code-block:: bash

   mkdir your_directory
   
However make sure to include pic.py, init_problem.py, the shock config file and any analytical scripts you wish to use in this directory.

Configuring the simulation
===================
To configure the simulation, you can edit the shock.ini file given. The following settings are given:

- [io]
   - outdir: Defines the output directory name
   - interval: No. of steps between output files written
   - full_interval: No. of steps between complete snapshot written to file (-1 = disabled)
   - restart: No. of steps between restart files written
   - stride: Output reduce factor
   - laprestart: When to restart (-1 = no restart, 0 = automatic, X lap to restart)

- [simulation]
   - cfl: The simulation time step in units of cfl
   - Nt: The number of laps in the simulation
   
- [problem]
   - Nspecies: Number of species
   - delgam: Temperature of plasma in units of mc^2 in bulk flow frame
   - temp_ratio: Ratio of ion temperature to electron temperature
   - gamma: Lorentz factor of bulk flow
   - me: Electron mass-to-charge
   - mi: ion mass-to-charge
   - sigma: Magnetization
   - npasses: Number of current filter passes
   - gammarad: Lorentz factor of radiation drag
   - radtemp: Temperature of radiation drag
   - bphi: B-field z angle
   - btheta: B-field x-y angle
   - wallgamma: X velocity of the left piston wall
   
- [grid]
   - Nx, Ny, Nz: No. of mesh grids in x,y,z direction
   - Nx/Ny/NzMesh: Size of mesh grids in x,y,z direction
   - dx/dy/dz: Size of spatial step
   
- [vmesh]
   - Nvx/y/z: Velocity mesh parameters

- [particles]
   - ppc: Particles per cell per species
   - c_omp: Speed of Light/Plasma Frequency
   - n_test_prtcls: Number of particles in simulation
   
Running the simulation
===================
To run a shock simulation on runko, use the following command:

.. code-block:: bash

   mpirun [-n no_of_cores] python3 pic.py --conf shock.ini

Using the Tools
===================
There are four tools provided for use in studying the results:

1. Particle Spectra
This script will output a spectra of the particles gamma over the time period of the simulation, called prtcl_spec.pdf. It will plot all laps by default however a lap can be designated. If all laps are plotted, it will also find the 10 most energetic particles at the end of the simulation and write their details to file in 10_prtcls.txt.

To run, use the command:

.. code-block:: bash

   python3 prtcl_spec.py --conf shock.ini [--lap lap_no]

2. Particle Path
This script will generate a file called prtcl_path.pdf which shows the 10 most energetic particles' position against time and gamma against time.

To run, use the command:

.. code-block:: bash

   python3 prtcl_path.py --conf shock.ini
  
3. Plot shock
This script generates a four-part graph showing:
- A density map
- A magnetic field map (Z-direction)
- A current density map (Z-direction)
- Plot of 1D density and magnetic energy density
This produces graphs for all laps generated, unless a lap is specified. Additionally, the paths of the 10 most energetic particles are shown and, if all laps are generated, the 1D density data is output for the shock velocity script.

To run, use the command:

.. code-block:: bash

   python3 plot_shock.py --conf shock.ini [--lap lap_no]

4. Shock velocity
This script will load in the 1D density data, and use it to find the compression ratio of the shock and subsequently the shock velocity based on the midpoint of the shock.
A graph of position against time will be shown, and the Beta-velocity value determined in the frame of the **downstream plasma**.

To run, use the command:

.. code-block:: bash

   python3 shock_RH.py --conf shock.ini
