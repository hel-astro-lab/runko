Tutorial
########

This tutorial is a walkthrough of different aspects constituting a runko simulation,
by going through pieces from turbulence example project at `projects/pic2-turbulence/pic2-decay.py`,
which is a self-contained version of decay setup from `projects/pic-turbulence/pic2.py`.
This tutorial focuses techical aspects and not on physics.

.. role:: python(code)
   :language: python


Configuration
=============

.. codeblock:: python
   :lineos:

   import runko

   config = runko.Configuration(None)

   config.Nt = 200 # Number of time steps.
   # How many tiles:
   config.Nx = 1
   config.Ny = 1
   config.Nz = 1
   # Size of emf grid in each tile:
   config.NxMesh = 80
   config.NyMesh = 80
   config.NzMesh = 80
   # ...


After importing `runko` we create a config object with `None` parameter.
This indicates that we want a empty config. `None` could be replaced by a path
to `.ini` file, but in order to keep the example self-contained,
we set the configuration parameters inline.

.. note::

   TODO: document configuration parameters


Logging
=======

.. codeblock:: python
   :lineos:

   logger = runko.runko_logger()
   # ...
   logger.info(f"q0: {config.q0}")
   logger.info(f"q1: {config.q1}")
   logger.info(f"m0: {config.m0}")
   logger.info(f"m1: {config.m1}")


Here we create a runko logger object and use its `info` method.
Notice that this will be executed on each MPI rank.
However, :python:`runko.runko_logger` is by default configured to write these to stdout
only from the main rank.

Many runko facilities use logging to inform what they are doing.
In addition to `info` method there is also `debug` logging method.
Output from it is by default disabled but can be enbaled by:

.. codeblock:: python
   :lineos:

    if runko.on_main_rank():
        import logging
        logger.setLevel(logging.DEBUG)


.. note::

   TODO: Write whole section about how the logging works in more detail
   and how it integrates with the Python standard library logging facilities.


Simulation initialization
=========================

Runko simulation is made up from 3D grid of tiles.
Each rank MPI rank owns some set of the tiles which are called local tiles.
How the tiles are distributed is configurable using `tile_partitioning` configuration parameter.

Communication between MPI ranks is done through so called virtual tiles.
They are the set of non-local tiles which are in moore neighborhood (i.e. includes diagonal neighbors)
of any local tile. Local tiles next to a virtual tiles are called boundary tiles.

Tiles are initialized through :python:`runo.TileGrid` object:

.. code:: python

   # TileGrid ctor:
   # - balances tiles based on conf (catepillar, hilbert)
   # - checks if restarts files are present for current config
   #   - if found initializes the tiles based on the restart files
   tile_grid = runko.TileGrid(config)

   if not tile_grid.initialized_from_restart_file():
       for idx in tile_grid.local_tile_indices():
           tile = runko.pic.Tile(idx, config)

           tile.batch_set_EBJ(zero_field, zero_field, zero_field,
                              Bx, By, Bz,
                              zero_field, zero_field, zero_field)

           # ppc = particles per cell (problem specific variable defined earlier)
           for _ in range(ppc):
               tile.batch_inject_to_cells(0, pgen0)
               tile.batch_inject_to_cells(1, pgen1)

           tile_grid.add_tile(tile, idx)


.. note::

   Restart files are note implemented in runko MVP.


We loop over indices corresponding to local tiles of this rank.
Then we construct a PIC tile and give it the grid index and configuration object.
For PIC tile we need to initialize electric field E, magnetic field B, current J and particles.
There are couple different ways to initialize them, but the most performant way is through
`tile.batch_set_EBJ` and `tile.batch_inject_to_cells`. These are explained below.
Lastly the initialized tile is added to the tile grid with `tile_grid.add_tile`
at specified tile index.

.. note::

   Runko MVP currently only supports a special case where each tile is the same,
   either `runko.pic.tile` or `runko.emf.tile`.
   Later there should be a way to have special kinds of tiles which can be used to implement
   e.g. boundary conditions other than periodic boundary condition.


Initializing fields
-------------------

`tile.batch_set_EBJ` takes nine parameters as input, one for each field component:
Ex, Ey, Ez, Bx, By, Bz, Jx, Jy and Jz.
The parameters are expected to be functions taking three parameters.
They are each invoked with three Numpy 3D `ndarray`s corresponding to grid coordinates x, y and z.
The functions should return a `ndarray` which is used to set the field values.

.. code:: python

   import numpy as np

   Bz = lambda x, y, z: np.full_like(x, binit)   # set each value to binit
   zero_field = lambda x, y, z: np.zeros_like(x) # set each value to zero

   def some_field(x, y, z):
       return np.sqrt(x**2 + y**2 + z**2) # set each value to distance from (0, 0, 0)


There is also simpler `tile.set_EBJ` which takes three functions, one for each field.
These functions are called with three floats and they should return the vector value
of the field at corresponding location as a tuple.


.. code:: python

   zero_field = lambda x, y, z: (0, 0, 0)
   tile.set_EBJ(zero_field, zero_field, zero_field)


Reason why `tile.batch_set_EBJ` takes nine parameters instead of three like `tile.set_EBJ`
is that the fields are defined on a Yee lattice, which means that different components
of a specific field are defined on slightly different locations,
so coordinates passed Bx, By and Bz functions are not exactly the same.
`tile.set_EBJ` goes around this by discarding the compontents which are evaluated at wrong locations.
With `tile.batch_set_EBJ` we don't want these inefficiencies.


Initializing particles
----------------------

`tile.batch_inject_to_cells` takes integer and a function as parameters.
The integer specifies which type of particle we are injecting.
Their mass and charge are defined by configuration parameters `mN` and `qN`,
where `N` is the particle type integer.

The given function is invoked with three 1D Numpy `ndarray`s
corresponding to all cell locations on a tile (cell's corner with smallest coordinates)
and it should return `runko.ParticleStateBatch`.

.. code:: python

    rng = np.random_default_rng(seed=42)

    # ...

    def pgen0(x, y, z):
        N = len(x)

        dx = rng.random(N)
        dy = rng.random(N)
        dz = rng.random(N)

        # Particles 1 are going on top of particles 0,
        # so these positions has to be saved such that pgen1 can get them.
        pgen0.pos = x + dx, y + dy, z + dz

        vel = runko.sample_boosted_juttner_synge(N, delgam0, beta=0, gen=rng)
        return runko.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    def pgen1(x, y, z):
        vel = runko.sample_boosted_juttner_synge(len(x), delgam1, beta=0, gen=rng)
        return runko.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    # ...

    tile.batch_inject_to_cells(0, pgen0)
    tile.batch_inject_to_cells(1, pgen1)


Runko uses units s.t. grid cells are 1 unit wide in each direction.
This allows to generate uniformly distributed positions inside cells with
:python:`x + dx, y + dy, z + dz`.

Velocities are generated according to boosted Jüttner-Synge distribution.
:python:`runko.sample_boosted_juttner_synge(N, ..., gen=rng)` generates `N`
3D velocities and uses the given Numpy random number generator.
It returns tuple of three 1D `ndarray`s.

Now `runko.ParticleStateBatch` can be constructed from position and velocity
tuples as shown in the code. For turbulence setup we want to generate type 1 particles
on top of type 0 particles. In `pgen0` we store the positions to a variable `pgen0.pos`
which can be accessed outside of the function and specifically from `pgen1`.

There exists simpler but slower ways to inject particles.
There is `tile.inject_to_each_cell` which differs from `tile.batch_inject_to_each_cell`
by calling the given function with three floats.
It also has to return `runko.ParticleState` which is like `runko.ParticleStateBatch`
but instead of position and velocity being tuples of `ndarray`s,
they are just tuples of three floats.
There is also `tile.inject` which takes integer particle type and a list of `runko.ParticleState`s.


Simulation
==========

After constructing the local tiles we can initialize the simulation with:

.. code:: python

   simulation = tile_grid.configure_simulation(config)

   def sync_EB(tile, comm, io):
       EB = (runko.comm_mode.emf_E, runko.comm_mode.emf_B)
       comm.virtual_tile_sync(*EB)
       comm.pairwise_moore(*EB)

       # Same as:
       # comm.virtual_tile_sync(runko.comm_mode.emf_B)
       # comm.virtual_tile_sync(runko.comm_mode.emf_E)
       # comm.pairwise_moore(runko.comm_mode.emf_B)
       # comm.pairwise_moore(runko.comm_mode.emf_E)

   simulation.prelude(sync_EB)


`tile_grid.configure_simulation` gives a handle to the simulation.
It is `runko.Simulation` object, but users should never try to construct it by hand.

Before the actual main simulation loop we do a single prelude step,
in order to not have a special case in the main loop for the first step.
Prelude step is defined using function which takes three opaque parameters.
Methods of `tile` are executed on each local tile and are specific to the used tile.
Methods of `comm` correspond to different differend kinds of communications
and lastly methods of `io` correspond to writing output.

.. code:: python

   def pic_simulation_step(tile, comm, io):

       tile.push_half_b()
       comm.virtual_tile_sync(runko.comm_mode.emf_B)
       comm.pairwise_moore(runko.comm_mode.emf_B)

       tile.push_particles()
       comm.virtual_tile_sync(runko.comm_mode.pic_particle)
       comm.pairwise_moore(runko.comm_mode.pic_particle)

       if simulation.lap % 5 == 0:
           tile.sort_particles()

       tile.deposit_current()
       comm.virtual_tile_sync(runko.comm_mode.emf_J)
       comm.pairwise_moore(runko.comm_mode.emf_J_exchange)

       comm.virtual_tile_sync(runko.comm_mode.emf_J)
       comm.pairwise_moore(runko.comm_mode.emf_J)
       tile.filter_current()
       comm.virtual_tile_sync(runko.comm_mode.emf_J)
       comm.pairwise_moore(runko.comm_mode.emf_J)
       tile.filter_current()
       tile.filter_current()

       tile.push_half_b()
       comm.virtual_tile_sync(runko.comm_mode.emf_B)
       comm.pairwise_moore(runko.comm_mode.emf_B)

       tile.push_e()
       tile.subtract_J_from_E()
       comm.virtual_tile_sync(runko.comm_mode.emf_E)
       comm.pairwise_moore(runko.comm_mode.emf_E)

       if simulation.lap % 20 == 0:
           io.emf_snapshot()

       if simulation.lap % 10 == 0:
           simulation.log_timer_statistics()


    simulation.for_each_lap(pic_simulation_step)
    simulation.log_timer_statistics()


The main simulation loop is executed with `simulation.for_each_lap`.
It will execute the given lap function while `simulation.lap` is less than config parameter `Nt`
(there is also `simulation.for_one_lap`).

Simulation automatically measures execution time of each step/method of the loop.
This information can be logged using `simulation.log_timer_statistics`.
For timer purposes each step are named based on the method name.
If there is multiple calls to the same method,
a running numbering is appended at the end.
For any method there is a possibility of explicitly naming them with `name` kwarg
(e.g. :python:`tile.pairwise_moore(runko.comm_mode.emf_B, name="foobar")`).


Communication
-------------

As an example let's look at only the beginning of the main loop
and see step by step what different communications do and why they are needed.
But first we have to understand what kind of halo regions runko uses.

Pic tile is refinment of emf tile which contains the underlying Yee lattice.
This lattice has to extend little bit further than actually belong to the tile,
in order to calculate derivatives at the boundaries
and to handle particles which are just on the borders of tiles.
This extended region is called the halo region and
runko makes practical choice of having three cells wide halo region to each direction.

.. code:: python

   def pic_simulation_step(tile, comm, io):

       tile.push_half_b()
       comm.virtual_tile_sync(runko.comm_mode.emf_B)
       comm.pairwise_moore(runko.comm_mode.emf_B)

       tile.push_particles()
       comm.virtual_tile_sync(runko.comm_mode.pic_particle)
       comm.pairwise_moore(runko.comm_mode.pic_particle)


All the `tile.*` methods that update the fields actually only update the non-halo region [1]_.
Therefor, after `tile.push_half_b` each tile has an outdated `B` in the halo region,
which prevents us executing `tile.push_particles` immediatly.

.. [1] With some exceptions. However, it is safest to always assume this.

In order to update it with the data from neighboring tile,
we can use `comm.pairwise_moore(runko.comm_mode.emf_B)`.
Pairwise refers to the fact that it consists many communications between tile pairs.
Moore refers to that the communication is done between local tiles
and tiles in their Moore neighborhood.
There is need for different kinds of pairwise Moore communications
but overall the shape of the communication stays the same:

- `runko.comm_mode.emf_{E,B,J}` modes just update halo regions with data from neighboring tile
   corresponding to the halo regions.
- `runko.comm_mode.pic_particle` transfers particles inside the halo to the corresponding
  tiles.
- `runko.comm_mode.emf_J_exchange` adds deposited current from boundary region of neighboring tile
  to the corresponding non-halo region of the "operated tile".
  This is needed as `tile.deposit_current` might generate current to halo-regions
  which has to be transfered to corresponding non-halo region.

Before we can actually do pairwise moore communication we have to sync the virtual tiles.
If we don't then boundary tiles will do pairwise communication with virtual tile
which is out of date with the corresponding tile on some other rank.
`comm.virtual_tile_sync` will update virtual tiles based on their corresponding "real" tiles.
Virtual tile sync with `runko.comm_mode.emf_{E,B,J}` will sync corresponding field data
and with `runko.comm_mode.particle` the particle data.

.. note::

   Here is a potential optimization.
   Instead of sending whole field/particle data in virtual tile sync,
   we could only send the actual data that is needed in some pairwise moore communication.
