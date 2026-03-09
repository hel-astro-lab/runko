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

    conf = runko.Configuration(args.conf)

    # --------------------------------------------------
    # Problem specific configuration

    ppc = conf.ppc
    sigma = conf.sigma
    c_omp = conf.c_omp
    theta = conf.theta                # dimensionless temperature kT/(mc^2)
    theta_ratio = conf.theta_ratio    # T_ion / T_electron
    upstream_gamma = conf.upstream_gamma
    n_filter_passes = conf.n_filter_passes or 3

    # Bulk flow (upstream_gamma is always a Lorentz factor >= 1)
    beta = np.sqrt(1.0 - 1.0 / upstream_gamma**2)

    # Charges and masses
    oppc = 2 * ppc
    omp = conf.cfl / c_omp
    conf.q0 = -(omp**2 * upstream_gamma) / (0.5 * oppc * (1.0 + conf.m0 / conf.m1))
    conf.q1 = abs(conf.q0)

    m0 = conf.m0 * abs(conf.q0)
    m1 = conf.m1 * abs(conf.q1)

    # Temperatures per species
    theta0 = theta                  # electron
    theta1 = theta_ratio * theta0   # ion/positron

    # B-field from sigma
    binit = np.sqrt(upstream_gamma * oppc * 0.5 * conf.cfl**2 * m0 * (1.0 + m0/m1) * sigma)

    # --------------------------------------------------
    # Field initialization
    Z  = lambda x, y, z: np.zeros_like(x)
    Bx = lambda x, y, z: np.full_like(x,       binit*conf.bx_proj)
    By = lambda x, y, z: np.full_like(x,       binit*conf.by_proj)
    Bz = lambda x, y, z: np.full_like(x,       binit*conf.bz_proj)
    Ey = lambda x, y, z: np.full_like(x, -beta*binit*conf.bz_proj)
    Ez = lambda x, y, z: np.full_like(x, +beta*binit*conf.by_proj)

    # --------------------------------------------------
    # Print setup summary
    if runko.on_main_rank():
        from runko.auto_outdir import resolve_outdir
        W = 21
        logger.info(f"{'--- [io] ---':}")
        logger.info(f"  {'outdir':<{W}}= {conf.outdir}")
        logger.info(f"  {'resolved outdir':<{W}}= {resolve_outdir(conf)}")
        logger.info(f"  {'prefix / postfix':<{W}}= {conf.prefix} / {conf.postfix}")

        logger.info(f"{'--- [grid] ---':}")
        logger.info(f"  {'tiles':<{W}}= {conf.Nx} x {conf.Ny} x {conf.Nz}")
        logger.info(f"  {'mesh per tile':<{W}}= {conf.NxMesh} x {conf.NyMesh} x {conf.NzMesh}")
        logger.info(f"  {'full grid':<{W}}= {conf.Nx*conf.NxMesh} x {conf.Ny*conf.NyMesh} x {conf.Nz*conf.NzMesh}")
        logger.info(f"  {'grid in c/wp':<{W}}= {conf.Nx*conf.NxMesh/c_omp:.1f} x {conf.Ny*conf.NyMesh/c_omp:.1f} x {conf.Nz*conf.NzMesh/c_omp:.1f}")

        logger.info(f"{'--- [simulation] ---':}")
        logger.info(f"  {'Nt':<{W}}= {conf.Nt}")
        logger.info(f"  {'cfl':<{W}}= {conf.cfl}")
        logger.info(f"  {'max plasma time':<{W}}= {conf.Nt * c_omp:.1f} wp^-1")

        logger.info(f"{'--- [particles] ---':}")
        logger.info(f"  {'ppc':<{W}}= {ppc}")
        logger.info(f"  {'m0 / m1':<{W}}= {conf.m0} / {conf.m1}")
        logger.info(f"  {'q0 / q1':<{W}}= {conf.q0:.6g} / {conf.q1:.6g}")
        logger.info(f"  {'theta_e / theta_i':<{W}}= {theta0:.6g} / {theta1:.6g}")

        logger.info(f"{'--- [problem] ---':}")
        logger.info(f"  {'sigma':<{W}}= {sigma}")
        logger.info(f"  {'c/wp':<{W}}= {c_omp}")
        logger.info(f"  {'upstream gamma':<{W}}= {upstream_gamma}")
        logger.info(f"  {'upstream beta':<{W}}= {beta:.6g}")
        logger.info(f"  {'B_init':<{W}}= {binit:.6g}")
        logger.info(f"  {'B direction':<{W}}= ({conf.bx_proj}, {conf.by_proj}, {conf.bz_proj})")
        logger.info(f"  {'n_filter_passes':<{W}}= {n_filter_passes}")
        logger.info(f"  {'v_A':<{W}}= {np.sqrt(sigma / (1.0 + sigma)):.6g}")
        logger.info(f"  {'beta_ion / beta_e':<{W}}= {2.0*theta1 / (sigma*(m1/m0+1.0)/(m1/m0)):.6g} / {2.0*theta0 / (sigma*(1.0/m0+1.0)):.6g}")

        logger.info(f"{'--- [algorithms] ---':}")
        logger.info(f"  {'field_propagator':<{W}}= {conf.field_propagator}")
        logger.info(f"  {'particle_pusher':<{W}}= {conf.particle_pusher}")
        logger.info(f"  {'field_interpolator':<{W}}= {conf.field_interpolator}")
        logger.info(f"  {'current_depositer':<{W}}= {conf.current_depositer}")
        logger.info(f"  {'current_filter':<{W}}= {conf.current_filter}")
        logger.info(f"  {'tile_partitioning':<{W}}= {conf.tile_partitioning}")

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

    tile_grid = runko.TileGrid(conf)

    if not tile_grid.initialized_from_restart_file():
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, conf)
            tile.batch_set_EBJ(Z, Ey, Ez, Bx, By, Bz, Z, Z, Z)
            tile.register_reflector_wall(wall)
            for _ in range(ppc):
                tile.batch_inject_to_cells(0, pgen0)
                tile.batch_inject_to_cells(1, pgen1)
            tile_grid.add_tile(tile, idx)

    # --------------------------------------------------
    # confure and start simulation

    simulation = tile_grid.configure_simulation(conf)
    simulation.verbose_laps = False

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
        x.io_ram_usage()

        if simulation.lap % 20 == 0:
            x.io_emf_snapshot()
            simulation.log_timer_statistics()

    simulation.for_each_lap(pic_simulation_step)
    simulation.log_timer_statistics()
