import runko

if __name__ == "__main__":

    conf = runko.Configuration(config_file="./turb_small.ini")
    logger = runko.Logger(conf)

    logger.start_timer("total")
    logger.start_timer("init")

    logger.log("Initializing simulation...")

    # TileGrid ctor:
    # - balances tiles based on conf (catepillar, hilbert)
    # - checks if restarts files are present for current config
    #   - if found initializes the tiles based on the restart files
    tile_grid = runko.TileGrid(conf)

    if not grid.initialized_from_restart_file():
        logger.log("Initializing simulation...")

        # There could be preconfigured cases which are initialized in pic tile ctor.
        # One option is to have user defined initial fields, as show here.

        E = lambda x, y, z: (0, 0, 0)
        B = runko.some_initial_conditions(conf)
        J = lambda x, y, z: (0, 0, 0)

        for idx in tile_grid.local_tile_indices():
            tile = runko.pypic.Tile(conf) # Tile can deduce its coordinate extents from conf.
            tile.set_fields(E, B, J)
            tile_grid.add_tile(idx, tile)

    # Initializes the simulation and returns a handle to it:
    # - analyze and sync boundaries between mpi tasks/ranks/localities
    # - loads virtual tiles (is there benefit of doing this explicitly?)
    # - based on conf loads (or give them explicitly as arguments):
    #   - field propagator
    #   - particle pusher
    #   - field interpolator
    #   - current depositor
    #   - current filter
    #   - antenna
    #   - all different writers
    simulation = tile_grid.configure_simulation(conf)

    logger.stop_timer("init")
    logger.log("Initialization finised.")

    # Could be runko.pic_simulation_step
    def pic_simulation_step(sol, communicate, writers):
        communicate("b", "e", name="mpi_b0_e0") # Implicitly updates boundaries.

        sol.field_propagator("push_half_b", nhood="local", name="push_half_b1")
        communicate("b", name="mpi_b1")

        solvers.particle_pusher(particles="all")  # Uses configured interpolator.

        solvers.field_propagator("push_half_b", nhood="local", name="push_half_b2")
        communicate("b", name="mpi_b2")

        solvers.field_propagator("push_e", nhood="local")

        # Could this be below particle pusher?
        communicate("exchange_particles")

        solvers.current_calculator(particles="all", nhood="local") # Automatically adds to E?
        communicate("exchange_currents")

        solvers.current_filter(nhood="local")
        solvers.antenna(nhood="local")

        # io through writers...

    simulation.for_each_lap(pic_simulation_step)

    logger.stop("total")
