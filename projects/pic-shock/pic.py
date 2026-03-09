"""
Particle-in-cell simulation of a collisionless shock
"""

import runko
import numpy as np
import argparse
import logging


if __name__ == "__main__":
    rng = np.random.default_rng(seed=42)

    logger = runko.runko_logger()

    if runko.on_main_rank():
        pass
        # logger.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(description="PIC collisionless shock simulation")
    parser.add_argument("--conf", type=str, default="shock_3d.ini",
                        help="Path to .ini configuration file")
    args = parser.parse_args()

    config = runko.Configuration(args.conf)

    # --------------------------------------------------
    # Problem specific configuration

    ppc = config.ppc
    sigma = config.sigma
    c_omp = config.c_omp
    theta = config.theta                # dimensionless temperature kT/(mc^2)
    theta_ratio = config.theta_ratio    # T_ion / T_electron
    upstream_gamma = config.upstream_gamma
    n_filter_passes = config.n_filter_passes or 3

    # Bulk flow (upstream_gamma is always a Lorentz factor >= 1)
    beta = np.sqrt(1.0 - 1.0 / upstream_gamma**2)

    # Charges and masses
    oppc = 2 * ppc
    omp = config.cfl / c_omp
    config.q0 = -(omp**2 * upstream_gamma) / (0.5 * oppc * (1.0 + config.m0 / config.m1))
    config.q1 = abs(config.q0)

    m0 = config.m0 * abs(config.q0)
    m1 = config.m1 * abs(config.q1)

    # Temperatures per species
    theta0 = theta                  # electron
    theta1 = theta_ratio * theta0   # ion/positron

    # B-field from sigma
    binit = np.sqrt(upstream_gamma * oppc * 0.5 * config.cfl**2 * m0 * (1.0 + m0/m1) * sigma)

    logger.info(f"Positron thermal spread: {theta1}")
    logger.info(f"Electron thermal spread: {theta0}")
    logger.info(f"Alfven vel: {np.sqrt(sigma / (1. + sigma))}")
    ion_beta = 2. * theta1 / (sigma * (m1/m0 + 1.) / (m1/m0))
    logger.info(f"Ion beta: {ion_beta}")
    logger.info(f"Electron beta: {2.*theta0/(sigma*(1/m0+1.))}")
    logger.info(f"sigma: {sigma}")
    logger.info(f"B_guide: {binit}")
    logger.info(f"q0: {config.q0}")
    logger.info(f"q1: {config.q1}")
    logger.info(f"m0: {config.m0}")
    logger.info(f"m1: {config.m1}")

    # --------------------------------------------------
    # Field initialization

    bx_proj = config.bx_proj or 0.0
    by_proj = config.by_proj or 0.0
    bz_proj = config.bz_proj or 1.0

    Bx0, By0, Bz0 = binit * bx_proj, binit * by_proj, binit * bz_proj
    Ey0, Ez0 = -beta * Bz0, +beta * By0

    Z  = lambda x, y, z: np.zeros_like(x)
    Bx = lambda x, y, z: np.full_like(x, Bx0)
    By = lambda x, y, z: np.full_like(x, By0)
    Bz = lambda x, y, z: np.full_like(x, Bz0)
    Ey = lambda x, y, z: np.full_like(x, Ey0)
    Ez = lambda x, y, z: np.full_like(x, Ez0)

    logger.info(f"Bx0: {Bx0}, By0: {By0}, Bz0: {Bz0}")
    logger.info(f"Ey0: {Ey0}, Ez0: {Ez0}")

    # --------------------------------------------------
    # Reflector wall (stationary for now)

    walloc = 15.0  # cell location of wall (leave >10 cells for BCs)
    wall = runko.pic.threeD.reflector_wall(walloc=walloc)
    logger.info(f"Reflector wall at x={walloc}")

    # --------------------------------------------------
    # Particle generators

    def pgen0(x, y, z):
        N = len(x)
        dx, dy, dz = rng.random(N), rng.random(N), rng.random(N)
        # Save positions so pgen1 places positrons on top of electrons
        pgen0.pos = x + dx, y + dy, z + dz
        vel = runko.sample_boosted_juttner_synge(
            N, theta0, Gamma=upstream_gamma, direction="-x", gen=rng)
        return runko.pic.threeD.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    def pgen1(x, y, z):
        vel = runko.sample_boosted_juttner_synge(
            len(x), theta1, Gamma=upstream_gamma, direction="-x", gen=rng)
        return runko.pic.threeD.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    # --------------------------------------------------
    # Grid and tile setup

    tile_grid = runko.TileGrid(config)

    if not tile_grid.initialized_from_restart_file():
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            tile.batch_set_EBJ(Z, Ey, Ez, Bx, By, Bz, Z, Z, Z)
            tile.register_reflector_wall(wall)
            for _ in range(ppc):
                tile.batch_inject_to_cells(0, pgen0)
                tile.batch_inject_to_cells(1, pgen1)
            tile_grid.add_tile(tile, idx)

    # --------------------------------------------------
    # Configure and start simulation

    simulation = tile_grid.configure_simulation(config)

    def sync_EB(x):
        EB = (runko.tools.comm_mode.emf_E, runko.tools.comm_mode.emf_B)
        x.comm_external(*EB)
        x.comm_local(*EB)

    simulation.prelude(sync_EB)

    # --------------------------------------------------
    # Main PIC loop

    def pic_simulation_step(x):

        # --- half B push + wall BC ---
        x.grid_push_half_b()
        x.grid_reflector_wall_field_bc()
        x.comm_external(runko.tools.comm_mode.emf_B)
        x.comm_local(runko.tools.comm_mode.emf_B)

        # --- particle push + reflect + communicate ---
        x.prtcl_push()
        x.prtcl_reflect_particles()
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

        if simulation.lap % 5 == 0:
            x.prtcl_sort()

        # --- current deposit + communicate ---
        x.prtcl_deposit_current()
        x.comm_external(runko.tools.comm_mode.emf_J)
        x.comm_local(runko.tools.comm_mode.emf_J_exchange)
        x.comm_external(runko.tools.comm_mode.emf_J)
        x.comm_local(runko.tools.comm_mode.emf_J)

        # --- current filter ---
        # Optimal: communicate every halo_size (=3) passes.
        # The separable binomial2 filter has a radius-1 stencil and writes to
        # [1, Nx-1), leaving boundary cells as zero. Each pass corrupts 1 more
        # halo cell inward from each side. With halo_size=3, up to 3 consecutive
        # passes keep the interior valid. After that, halos must be refreshed.
        for i in range(n_filter_passes):
            if i > 0 and i % 3 == 0:
                x.comm_external(runko.tools.comm_mode.emf_J)
                x.comm_local(runko.tools.comm_mode.emf_J)
            x.grid_filter_current()

        # --- second half B push + wall BC ---
        x.grid_push_half_b()
        x.grid_reflector_wall_field_bc()
        x.comm_external(runko.tools.comm_mode.emf_B)
        x.comm_local(runko.tools.comm_mode.emf_B)

        # --- E push + wall BC + add current ---
        x.grid_push_e()
        x.grid_reflector_wall_field_bc()
        x.grid_add_current()
        x.grid_reflector_wall_field_bc()
        x.comm_external(runko.tools.comm_mode.emf_E)
        x.comm_local(runko.tools.comm_mode.emf_E)

        # --- advance wall position (for moving wall; no-op for stationary) ---
        x.prtcl_advance_reflector_walls()

        # --- IO ---
        x.io_average_kinetic_energy()
        x.io_average_B_energy_density()
        x.io_average_E_energy_density()

        if simulation.lap % 20 == 0:
            x.io_emf_snapshot()
            simulation.log_timer_statistics()

    simulation.for_each_lap(pic_simulation_step)
    simulation.log_timer_statistics()
