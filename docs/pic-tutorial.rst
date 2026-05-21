PIC tutorial
############

This tutorial walks through the different pieces that make up a runko PIC simulation,
using the PIC turbulence example project at ``projects/pic-turbulence/pic.py`` as the running example.
It focuses on technical aspects rather than on the underlying physics.

.. role:: python(code)
   :language: python


Configuration
=============

.. code:: python

   import runko

   config = runko.Configuration(None)

   config.Nt = 200 # Number of time steps.
   # How many tiles:
   config.Nx = 4
   config.Ny = 4
   config.Nz = 4
   # Size of emf grid in each tile:
   config.NxMesh = 20
   config.NyMesh = 20
   config.NzMesh = 20
   # ...


After importing :python:`runko` we create a config object with a :python:`None` parameter.
This indicates that we want an empty config. :python:`None` could be replaced by a path
to an ``.ini`` file, but to keep the example self-contained,
we set the configuration parameters inline.

.. note::

   TODO: document configuration parameters


Logging
=======

.. code:: python

   logger = runko.runko_logger()
   # ...
   logger.info(f"q0: {config.q0}")
   logger.info(f"q1: {config.q1}")
   logger.info(f"m0: {config.m0}")
   logger.info(f"m1: {config.m1}")


Here we create a runko logger object and use its :python:`info` method.
Note that this call runs on every MPI rank;
however, :python:`runko.runko_logger` is configured by default to write to stdout
only from the main rank.

Many runko facilities use logging to report what they are doing.
In addition to the :python:`info` method, there is also a :python:`debug` logging method.
Its output is disabled by default but can be enabled with:

.. code:: python

    if runko.on_main_rank():
        import logging
        logger.setLevel(logging.DEBUG)


.. note::

   TODO: Write whole section about how the logging works in more detail
   and how it integrates with the Python standard library logging facilities.


Simulation initialization
=========================

A runko simulation is made up of a 3D grid of tiles.
Each MPI rank owns some subset of the tiles, which are called local tiles.
How the tiles are distributed is configurable via the :python:`tile_partitioning` configuration parameter.

Communication between MPI ranks goes through so-called virtual tiles.
These are the non-local tiles that lie in the Moore neighborhood (i.e. including diagonal neighbors)
of any local tile. Local tiles adjacent to virtual tiles are called boundary tiles.

Tiles are initialized through a :python:`runko.TileGrid` object:

.. code:: python

   tile_grid = runko.TileGrid(config)

   if not tile_grid.initialized_from_restart_file():
       for idx in tile_grid.local_tile_indices():
           tile = runko.pic.threeD.Tile(idx, config)
           tile.batch_set_EBJ(zero_field, zero_field, zero_field,
                              Bx, By, Bz,
                              zero_field, zero_field, zero_field)

           # ppc = particles per cell (problem specific variable defined earlier)
           for _ in range(ppc):
                tile.batch_inject_to_cells(0, pgen0)
                tile.batch_inject_to_cells(1, pgen1)
           tile_grid.add_tile(tile, idx)


.. note::

   Restart files are not yet implemented.


We loop over the indices corresponding to local tiles of this rank.
For each one we construct a PIC tile, passing the grid index and the configuration object.
A PIC tile needs initial values for the electric field E, magnetic field B, current J, and particles.
There are a couple of different ways to initialize them, but the most performant route is via
:python:`tile.batch_set_EBJ` and :python:`tile.batch_inject_to_cells`. Both are explained below.
Finally, the initialized tile is added to the tile grid with :python:`tile_grid.add_tile`
at the specified tile index.

.. note::

   Runko currently only supports the special case where every tile has the same type,
   either :python:`runko.pic.tile` or :python:`runko.emf.tile`.
   In the future there should be a way to mix specialized tile types — for example,
   to implement boundary conditions other than periodic.


Initializing fields
-------------------

:python:`tile.batch_set_EBJ` takes nine parameters as input, one for each field component:
Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, and Jz.
Each parameter is expected to be a function taking three parameters.
They are each invoked with three NumPy 3D :python:`ndarray` s corresponding to the grid coordinates x, y, and z.
The functions should return an :python:`ndarray` whose values are used to set the field.

.. code:: python

   import numpy as np

   Bz = lambda x, y, z: np.full_like(x, binit)   # set each value to binit
   zero_field = lambda x, y, z: np.zeros_like(x) # set each value to zero

   def some_field(x, y, z):
       return np.sqrt(x**2 + y**2 + z**2) # set each value to distance from (0, 0, 0)


There is also a simpler :python:`tile.set_EBJ`, which takes three functions, one for each field.
These functions are called with three floats and should return the vector value
of the field at the corresponding location as a tuple.


.. code:: python

   zero_field = lambda x, y, z: (0, 0, 0)
   tile.set_EBJ(zero_field, zero_field, zero_field)


The reason :python:`tile.batch_set_EBJ` takes nine parameters
instead of three like :python:`tile.set_EBJ`
is that the fields are defined on a Yee lattice: different components
of a given field live at slightly different locations,
so the coordinates passed to the Bx, By, and Bz functions are not exactly the same.
:python:`tile.set_EBJ` works around this by discarding the components
that are evaluated at the wrong locations.
With :python:`tile.batch_set_EBJ` we avoid that overhead.


Initializing particles
----------------------

:python:`tile.batch_inject_to_cells` takes an integer and a function as parameters.
The integer specifies which species of particle is being injected.
A species' mass and charge are defined by the configuration parameters :python:`mN` and :python:`qN`,
where :python:`N` is the species integer.

The given function is invoked with three 1D NumPy :python:`ndarray` s
corresponding to every cell location on the tile (the corner of each cell with the smallest coordinates)
and it should return a :python:`runko.ParticleStateBatch`.

.. code:: python

   rng = np.random.default_rng(seed=42)

   # ...

   def pgen0(x, y, z):
       N = len(x)

       dx = rng.random(N)
       dy = rng.random(N)
       dz = rng.random(N)

       # Particles of species 1 will be placed on top of species 0,
       # so these positions have to be saved so that pgen1 can reuse them.
       pgen0.pos = x + dx, y + dy, z + dz

       vel = runko.sample_boosted_juttner_synge(N, delgam0, beta=0, gen=rng)
       return runko.ParticleStateBatch(pos=pgen0.pos, vel=vel)

   def pgen1(x, y, z):
       vel = runko.sample_boosted_juttner_synge(len(x), delgam1, beta=0, gen=rng)
       return runko.ParticleStateBatch(pos=pgen0.pos, vel=vel)

   # ...

    tile.batch_inject_to_cells(0, pgen0)
    tile.batch_inject_to_cells(1, pgen1)


Runko uses units such that grid cells are 1 unit wide in each direction.
This makes it possible to generate uniformly distributed positions inside cells with
:python:`x + dx, y + dy, z + dz`.

Velocities are drawn from a boosted Jüttner–Synge distribution.
:python:`runko.sample_boosted_juttner_synge(N, ..., gen=rng)` generates :python:`N`
3D velocities using the given NumPy random number generator.
It returns a tuple of three 1D :python:`ndarray` s.

A :python:`runko.ParticleStateBatch` can then be constructed from the position and velocity
tuples as shown in the code. For the turbulence setup we want to place species-1 particles
on top of species-0 particles;
to do this, :python:`pgen0` stores the positions in a function attribute :python:`pgen0.pos`
that can be accessed from the outside, and specifically from :python:`pgen1`.

There are simpler but slower ways to inject particles.
:python:`tile.inject_to_each_cell` differs from :python:`tile.batch_inject_to_cells`
in that it calls the given function with three floats instead of three arrays.
It has to return a :python:`runko.ParticleState`, which is like :python:`runko.ParticleStateBatch`
except that position and velocity are tuples of three floats rather than tuples of :python:`ndarray` s.
There is also :python:`tile.inject`, which takes an integer particle type
and a list of :python:`runko.ParticleState` objects.


Simulation
==========

After constructing the local tiles we can initialize the simulation with:

.. code:: python

   simulation = tile_grid.configure_simulation(config)

   def sync_EB(x):
       EB = (runko.tools.comm_mode.emf_E, runko.tools.comm_mode.emf_B)
       x.comm_external(*EB)
       x.comm_local(*EB)

   simulation.prelude(sync_EB)


:python:`tile_grid.configure_simulation` returns a handle to the simulation.
It is a :python:`runko.Simulation` object, but users should never try to construct it by hand.

Before the actual main simulation loop we run a single prelude step,
so that the main loop does not need a special case for the first iteration.
The prelude step is defined using a function that takes a single parameter.

.. code:: python

   def pic_simulation_step(x):

       x.grid_push_half_b()
       x.comm_external(runko.tools.comm_mode.emf_B)
       x.comm_local(runko.tools.comm_mode.emf_B)

       x.prtcl_push()
       x.prtcl_pack_outgoing()
       x.comm_external(runko.tools.comm_mode.pic_particle)
       x.comm_local(runko.tools.comm_mode.pic_particle)

       if simulation.lap % 5 == 0:
           x.prtcl_sort()

       x.prtcl_deposit_current()
       x.comm_external(runko.tools.comm_mode.emf_J)
       x.comm_local(runko.tools.comm_mode.emf_J_exchange)

       x.comm_external(runko.tools.comm_mode.emf_J)
       x.comm_local(runko.tools.comm_mode.emf_J)

       x.grid_filter_current()
       x.comm_external(runko.tools.comm_mode.emf_J)
       x.comm_local(runko.tools.comm_mode.emf_J)
       x.grid_filter_current()
       x.grid_filter_current()

       x.grid_push_half_b()
       x.comm_external(runko.tools.comm_mode.emf_B)
       x.comm_local(runko.tools.comm_mode.emf_B)

       x.grid_push_e()
       x.grid_add_current()
       x.comm_external(runko.tools.comm_mode.emf_E)
       x.comm_local(runko.tools.comm_mode.emf_E)

       x.io_average_kinetic_energy()

       if simulation.lap % 20 == 0:
           x.io_emf_snapshot()

       if simulation.lap % 20 == 0:
           simulation.log_timer_statistics()

   simulation.for_each_lap(pic_simulation_step)
   simulation.log_timer_statistics()


The main simulation loop is executed with :python:`simulation.for_each_lap`.
It runs the given lap function while :python:`simulation.lap` is less than
the config parameter :python:`Nt`
(there is also :python:`simulation.for_one_lap` for running just a single lap).

The simulation automatically measures the execution time of each step/method in the loop.
This information can be logged using :python:`simulation.log_timer_statistics`.
For timer purposes, each step is named after the method it calls.
If there are multiple calls to the same method,
a running index is appended to disambiguate them.
Any method can also be explicitly named with the :python:`name` kwarg
(e.g. :python:`tile.comm_external(runko.comm_mode.emf_B, name="foobar")`).


Communication
-------------

As an example, let us look at just the beginning of the main loop
and step through what the different communications do and why they are needed.
But first we have to understand what kind of halo regions runko uses.

A PIC tile is a refinement of an emf tile, which contains the underlying Yee lattice.
This lattice has to extend slightly beyond the tile's own cells,
both to compute derivatives at the boundaries
and to handle particles that sit right on the edge of a tile.
This extended region is called the halo region;
runko makes the practical choice of a three-cell-wide halo in each direction.

.. code:: python

   def pic_simulation_step(x):

       x.grid_push_half_b()
       x.comm_external(runko.tools.comm_mode.emf_B)
       x.comm_local(runko.tools.comm_mode.emf_B)

       x.prtcl_push()
       x.prtcl_pack_outgoing()
       x.comm_external(runko.tools.comm_mode.pic_particle)
       x.comm_local(runko.tools.comm_mode.pic_particle)


For example, :python:`x.push_half_b` only updates values in the non-halo region.
As a result, each tile is left with an outdated :python:`B` in its halo region,
which prevents us from executing :python:`x.prtcl_push` immediately.

To refresh the halo with data from neighboring tiles,
we use :python:`x.comm_local(runko.comm_mode.emf_B)`.
Several different kinds of local communication exist,
but the overall pattern stays the same:

- :python:`runko.comm_mode.emf_{E,B,J}` modes simply update halo regions with the corresponding
  data from neighboring tiles.
- :python:`runko.comm_mode.pic_particle` transfers particles inside the halo to the
  corresponding tiles.
- :python:`runko.comm_mode.emf_J_exchange` adds deposited current from the halo region of a neighboring tile
  to the corresponding non-halo region of the "operated tile".
  This is needed because :python:`tile.deposit_current` may generate current in the halo region,
  which then has to be transferred to the corresponding non-halo region.

Before we can do the local communication, we have to sync the virtual tiles.
Otherwise, boundary tiles will perform their local communication against a virtual tile
that is out of date with respect to the corresponding tile on another rank.
:python:`x.comm_external` updates virtual tiles based on their corresponding "real" tiles.
External communication with :python:`runko.comm_mode.emf_{E,B,J}` syncs the corresponding field data,
and with :python:`runko.comm_mode.particle` it syncs the particle data.
Before particle communication, outgoing particles have to be packed with :python:`prtcl_pack_outgoing`.


.. note::

   TODO: document I/O
