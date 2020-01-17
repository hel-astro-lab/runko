Tutorial: Collisionless shocks
##############################

In this tutorial we use Runko to run a small nonmagnetized collisionless shock simulation.

.. image:: https://cdn.jsdelivr.net/gh/natj/pb-utilities@master/movies/shock.gif


Running a collisionless shock simulation
========================================
You will find the related scripts to run the simulation and the analytical tools in

.. code-block:: bash

   /runko/projects/shocks


Configuring the simulation
==========================
To configure the simulation, you can edit the `shock_mini.ini` or `shock_small.ini` files given. The smaller of those (mini) should finish in less than a minute whereas the small takes ~20mins to finish.

The following settings are given:

- [io]
   - `outdir`: Defines the output directory name
   - `interval`: No. of steps between output files written
   - `restart`: No. of steps between restart files written
   - `laprestart`: Simulation lap to restart from (-1 = no restart, 0 = automatic, X lap to restart)

Additionally, more advanced options include

   - `stride`: Output image reduce factor
   - `full_interval`: No. of steps between full restart snapshot written to file (-1 = disabled)

- [simulation]
   - `cfl`: The simulation time step in units of grid space
   - `Nt`: The total number of laps in the simulation
   
- [problem]
   - `Nspecies`: Number of species
   - `delgam`: Temperature of plasma in units of mc^2 in bulk flow frame
   - `gamma`: Lorentz factor of bulk flow
   - `me`: electron mass-to-charge
   - `mi`: positron/ion mass-to-charge
   - `sigma`: Plasma magnetization parameter
   - `npasses`: Number of current filter passes
   
Additionally, more advanced options include

   - `temp_ratio`: Ratio of ion temperature to electron temperature
   - `bphi`: external B-field z angle
   - `btheta`: external B-field x-y angle
   - `wallgamma`: X velocity of the left piston wall (this makes the left wall ram into the plasma)

- [grid]
   - `Nx`, `Ny`, `Nz`: No. of mesh tiles in x,y,z direction
   - `NxMesh`/`NyMesh`/`NzMesh`: Size of tiles in x,y,z direction
   
Complete grid size is then `Ni*NiMesh`.

- [particles]
   - `ppc`: Particles per cell per species
   - `c_omp`: Skindepth resolution

Additionally, more advanced options include

   - `n_test_prtcls`: Number of test particles tracked in simulation
   


Running the simulation
======================
To run a shock simulation on Runko, use the following command:

.. code-block:: bash

   mpirun [-n no_of_cores] python3 pic.py --conf shock_mini.ini


Using the analysis tools
========================
There are four tools provided for use in studying the results:

1. Particle spectra
-------------------

This script will output a spectra of the particles Lorentz factor (gamma) over the time period of the simulation. It will plot all simulation laps by default, however, a lap can also be designated. If all laps are plotted, it will also find the 10 most energetic particles at the end of the simulation and write their details to file in `10_prtcls.txt`.

To run, use the command:

.. code-block:: bash

   python3 prtcl_spec.py --conf shock_mini.ini [--lap lap_no]

2. Particle paths
-----------------

This script will generate a file which shows the 10 most energetic particles' position and their Lorentz factor against time.

To run, use the command:

.. code-block:: bash

   python3 prtcl_path.py --conf shock_mini.ini
  
3. Plot shock
-------------

This script generates a four-part plot showing:

- A plasma density map
- An out-of-the-plane magnetic field map (Z-direction)
- An out-of-the-plane current density map (Z-direction)
- Plot of 1D density and magnetic energy density

This produces graphs for all laps generated, unless a simulation lap is specified. Additionally, the paths of the 10 most energetic particles are shown and, if all laps are generated, the 1D density data is output for the shock velocity script.

To run, use the command:

.. code-block:: bash

   python3 plot_shock.py --conf shock_mini.ini [--lap lap_no]

4. Shock Rankine-Hugoniot jump conditions
-----------------------------------------

This script will load in the 1D density data, and use it to find the compression ratio of the shock and subsequently the shock velocity based on the midpoint of the shock.
A plot of position against time will be shown, and the coordinate-velocity value determined in the frame of the **downstream plasma**.

To run, use the command:

.. code-block:: bash

   python3 shock_RH.py --conf shock_mini.ini
