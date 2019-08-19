Running the shock and tools
############

In this section it is explained how to run Runko and the given analytical tools

Using the script
===================
A python script is provided called Run.py, which can automatically run the shock and analytical tools, along with ffmpeg in succession. Details of each are provided below.
To run the script, use the command:
.. code-block:: bash

   python3 Run.py --conf shock.ini
   
It will ask if you wish to run the analytical tools and ffmpeg or just the shock, and how many cores you wish to use in the simulation.

Using the Tools
===================
There are four tools provided for use in studying the results:

1. Particle Spectra
This script will output a spectra of the particles gamma over the time period of the simulation, called prtcl_spec.pdf. It will plot all laps by default however a lap can be designated. If all laps are plotted, it will also find the 10 most energetic particles at the end of the simulation and write their details to file in 10_prtcls.txt.
To run, use the command:

.. code-block:: bash

   python3 prtcl_spec.py --conf shock.ini [--lap lap_no]

2. Particle Path
This script will generate a file called prtcl_path.png which shows the 10 most energetic particles' position against time and gamma against time.
To run, use the command:

.. code-block:: bash

   python3 prtcl_path.py --conf shock.ini
  
3. UniPlot
This script generates a four-part graph showing:
- A density map
- A magnetic field map (Z-direction)
- A current density map (Z-direction)
- Plot of 1D density and magnetic energy density
This produces graphs for all laps generated, unless a lap is specified. Additionally, the paths of the 10 most energetic particles are shown and, if all laps are generated, the 1D density data is output for the shock velocity script.
To run, use the command:

.. code-block:: bash

   python3 UniPlot.py --conf shock.ini [--lap lap_no]

4. Shock velocity
This script will load in the 1D density data, and use it to find the compression ratio of the shock and subsequently the shock velocity based on the midpoint of the shock.
A graph of position against time will be shown, and the Beta-velocity value determined in the frame of the **downstream plasma**.
To run, use the command:

.. code-block:: bash

   python3 RH_shk.py --conf shock.ini
