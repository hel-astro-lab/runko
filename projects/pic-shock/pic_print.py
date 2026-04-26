"""
Particle-in-cell simulation of a collisionless shock

DEBUG VERSION of pic.py: every call inside the per-lap step is preceded by
an MPI.COMM_WORLD.Barrier() and a rank-0 announce. The last printed
"[lap N] -> <name>" line names the call that segfaulted.

Instrumentation is gated by DEBUG_FROM_LAP (hard-coded below) so noisy
ramp-up laps don't flood the log; only laps >= DEBUG_FROM_LAP are traced.
"""

import runko
import numpy as np
import argparse
import logging

from mpi4py import MPI


if __name__ == "__main__":
    rng = np.random.default_rng(seed=42)

    logger = runko.runko_logger()

    if runko.on_main_rank():
        pass
        # logger.setLevel(logging.DEBUG)

    # ----------------------------------------------------------------
    # Debug instrumentation: barrier + rank-0 announce before each call.
    # Tune DEBUG_FROM_LAP to skip ramp-up; 0 = trace every lap.
    _world = MPI.COMM_WORLD
    DEBUG_FROM_LAP = 150

    def _sync_and_announce(name):
        if simulation.lap < DEBUG_FROM_LAP:
            return
        _world.Barrier()
        if runko.on_main_rank():
            logger.info(f"[lap {simulation.lap}] -> {name}")
    # ----------------------------------------------------------------

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
    output_interval = conf.output_interval or 20

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
    # Upstream field values
    Ex_up = 0.0
    Ey_up = -beta*binit*conf.bz_proj
    Ez_up = +beta*binit*conf.by_proj
    Bx_up = binit*conf.bx_proj
    By_up = binit*conf.by_proj
    Bz_up = binit*conf.bz_proj

    # Field initialization
    Z  = lambda x, y, z: np.zeros_like(x)
    Bx = lambda x, y, z: np.full_like(x, Bx_up)
    By = lambda x, y, z: np.full_like(x, By_up)
    Bz = lambda x, y, z: np.full_like(x, Bz_up)
    Ey = lambda x, y, z: np.full_like(x, Ey_up)
    Ez = lambda x, y, z: np.full_like(x, Ez_up)

    # --------------------------------------------------
    # Print setup summary
    if runko.on_main_rank():
        from runko.auto_outdir import resolve_outdir
        W = 21
        logger.info(f"{'--- [io] ---':}")
        logger.info(f"  {'outdir':<{W}}= {conf.outdir}")
        logger.info(f"  {'resolved outdir':<{W}}= {resolve_outdir(conf)}")
        logger.info(f"  {'prefix / postfix':<{W}}= {conf.prefix} / {conf.postfix}")
        logger.info(f"  {'output_interval':<{W}}= {output_interval}")
        logger.info(f"  {'spectra_nbins':<{W}}= {conf.spectra_nbins}")
        logger.info(f"  {'spectra_umin':<{W}}= {conf.spectra_umin}")
        logger.info(f"  {'spectra_umax':<{W}}= {conf.spectra_umax}")

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
        logger.info(f"  {'Lx crossing':<{W}}= {conf.Nx*conf.NxMesh/conf.cfl:.0f} laps")
        logger.info(f"  {'beta_ion / beta_e':<{W}}= {2.0*theta1 / (sigma*(m1/m0+1.0)/(m1/m0)):.6g} / {2.0*theta0 / (sigma*(1.0/m0+1.0)):.6g}")

        logger.info(f"{'--- [algorithms] ---':}")
        logger.info(f"  {'field_propagator':<{W}}= {conf.field_propagator}")
        logger.info(f"  {'particle_pusher':<{W}}= {conf.particle_pusher}")
        logger.info(f"  {'field_interpolator':<{W}}= {conf.field_interpolator}")
        logger.info(f"  {'current_depositer':<{W}}= {conf.current_depositer}")
        logger.info(f"  {'current_filter':<{W}}= {conf.current_filter}")
        logger.info(f"  {'tile_partitioning':<{W}}= {conf.tile_partitioning}")

        logger.info(f"--- [debug] ---")
        logger.info(f"  {'DEBUG_FROM_LAP':<{W}}= {DEBUG_FROM_LAP}")

    # --------------------------------------------------
    # Reflector wall and moving injector

    Lx = conf.Nx * conf.NxMesh
    walloc = 15.0  # cell location of wall (leave >10 cells for BCs)
    wall = runko.pic.threeD.reflector_wall(walloc=walloc)
    conducting_bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=walloc,
                                              E_components=0b110, B_components=0, J_components=0b111)
    upstream_bc = runko.emf.threeD.edge_bc(direction=0, side=1, position=Lx - 5,
                                            Ex=Ex_up, Ey=Ey_up, Ez=Ez_up,
                                            Bx=Bx_up, By=By_up, Bz=Bz_up,
                                            J_components=0b111)

    injloc0 = walloc + 10.0 * c_omp  # initial injection right edge
    n_inj = 50
    injector = runko.MovingInjector(
        injloc=injloc0, beta_inj=1.0, beta_flow=beta,
        cfl=conf.cfl, n_inj=n_inj, walloc=walloc, Lx=Lx)

    logger.info(f"Reflector wall at x={walloc}")
    logger.info(f"Moving injector at x={injloc0}, beta_inj=1.0, beta_flow={beta:.4g}, n_inj={n_inj}")

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

    # Particle container pre-allocation. Two optional knobs (either or both
    # may be set in [particles]; taking the tighter constraint):
    #   prealloc_factor            — multiplier on ppc (e.g. 4.0 for the
    #                                Rankine-Hugoniot downstream compression
    #                                ratio). 0 / unset = no pre-allocation.
    #   prealloc_memory_gb_per_rank — cap by per-rank memory budget
    #                                (32 B / particle, split across local
    #                                tiles and species).
    n_species = 2
    local_tiles_count = sum(1 for _ in tile_grid.local_tile_indices())
    cells_per_tile = conf.NxMesh * conf.NyMesh * conf.NzMesh
    prealloc_caps = []
    if conf.prealloc_factor is not None:
        prealloc_caps.append(int(conf.prealloc_factor * ppc * cells_per_tile))
    if conf.prealloc_memory_gb_per_rank is not None and local_tiles_count > 0:
        budget_bytes = conf.prealloc_memory_gb_per_rank * 1e9
        prealloc_caps.append(int(budget_bytes / local_tiles_count / n_species / 32))
    prealloc_n = min(prealloc_caps) if prealloc_caps else 0

    # Pass per-species prealloc size through the config so each Tile sizes its
    # particle containers at construction time with dead-filled slots; later
    # batch_inject_in_x_stripe calls write into those slots in place without
    # reallocation.
    conf.prealloc_per_species = prealloc_n

    if runko.on_main_rank() and prealloc_n > 0:
        logger.info(
            f"Pre-allocating {prealloc_n} dead particles per species per tile "
            f"at construction ({prealloc_n * n_species * 32 / 1e6:.1f} MB/tile, "
            f"{prealloc_n * n_species * 32 * local_tiles_count / 1e9:.2f} GB/rank)"
        )

    if True: # regular shock setup
        if not tile_grid.initialized_from_restart_file():
            for idx in tile_grid.local_tile_indices():
                tile = runko.pic.threeD.Tile(idx, conf)
                tile.batch_set_EBJ(Z, Ey, Ez, Bx, By, Bz, Z, Z, Z)
                tile.register_reflector_wall(wall)
                tile.register_edge_bc(conducting_bc)
                tile.register_edge_bc(upstream_bc)
                for _ in range(ppc):
                    tile.batch_inject_in_x_stripe(0, pgen0, walloc, injloc0)
                    tile.batch_inject_in_x_stripe(1, pgen1, walloc, injloc0)
                tile_grid.add_tile(tile, idx)

    else: # Cherenkov-measuring setup w/o reflectors or boundaries
        if not tile_grid.initialized_from_restart_file():
            for idx in tile_grid.local_tile_indices():
                tile = runko.pic.threeD.Tile(idx, conf)
                tile.batch_set_EBJ(Z, Ey, Ez, Bx, By, Bz, Z, Z, Z)
                # no reflectors or bcs
                for _ in range(ppc):
                    # inject everywhere
                    tile.batch_inject_to_cells(0, pgen0)
                    tile.batch_inject_to_cells(1, pgen1)
                tile_grid.add_tile(tile, idx)


    # --------------------------------------------------
    # Configure and start simulation

    simulation = tile_grid.configure_simulation(conf)
    simulation.verbose_laps = False # do not print lap info every time; only statistics

    def sync_EB(x):
        EB = (runko.tools.comm_mode.emf_E, runko.tools.comm_mode.emf_B)
        x.comm_external(*EB)
        x.comm_local(*EB)

    simulation.prelude(sync_EB)

    # --------------------------------------------------
    # Main PIC loop (instrumented)

    def pic_simulation_step(x):

        # --- half B push + B edge BCs ---
        _sync_and_announce("grid_push_half_b [1st]")
        x.grid_push_half_b()
        _sync_and_announce("grid_apply_edge_bcs(emf_B) [after 1st half B]")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_B)
        _sync_and_announce("comm_external(emf_B) [after 1st half B]")
        x.comm_external(runko.tools.comm_mode.emf_B)
        _sync_and_announce("comm_local(emf_B) [after 1st half B]")
        x.comm_local(runko.tools.comm_mode.emf_B)

        # --- particle push + reflect + communicate ---
        _sync_and_announce("prtcl_push")
        x.prtcl_push()
        _sync_and_announce("prtcl_reflect_particles")
        x.prtcl_reflect_particles()
        _sync_and_announce("prtcl_pack_outgoing")
        x.prtcl_pack_outgoing()
        _sync_and_announce("comm_external(pic_particle)")
        x.comm_external(runko.tools.comm_mode.pic_particle)
        _sync_and_announce("comm_local(pic_particle)")
        x.comm_local(runko.tools.comm_mode.pic_particle)

        if simulation.lap % 5 == 0:
            _sync_and_announce("prtcl_sort")
            x.prtcl_sort()

        # --- current deposit + communicate ---
        _sync_and_announce("prtcl_deposit_current")
        x.prtcl_deposit_current()
        _sync_and_announce("comm_external(emf_J) [post-deposit #1]")
        x.comm_external(runko.tools.comm_mode.emf_J)
        _sync_and_announce("comm_local(emf_J_exchange)")
        x.comm_local(runko.tools.comm_mode.emf_J_exchange)
        _sync_and_announce("comm_external(emf_J) [post-deposit #2]")
        x.comm_external(runko.tools.comm_mode.emf_J)
        _sync_and_announce("comm_local(emf_J) [post-deposit]")
        x.comm_local(runko.tools.comm_mode.emf_J)
        _sync_and_announce("grid_apply_edge_bcs(emf_J) [post-deposit]")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_J)

        # --- current filter ---
        for i in range(n_filter_passes):
            if i > 0 and i % 3 == 0:
                _sync_and_announce(f"comm_external(emf_J) [filter pass {i} refresh]")
                x.comm_external(runko.tools.comm_mode.emf_J)
                _sync_and_announce(f"comm_local(emf_J) [filter pass {i} refresh]")
                x.comm_local(runko.tools.comm_mode.emf_J)
            _sync_and_announce(f"grid_filter_current [pass {i}]")
            x.grid_filter_current()
        _sync_and_announce("grid_apply_edge_bcs(emf_J) [post-filter]")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_J)

        # --- second half B push + B edge BCs ---
        _sync_and_announce("grid_push_half_b [2nd]")
        x.grid_push_half_b()
        _sync_and_announce("grid_apply_edge_bcs(emf_B) [after 2nd half B]")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_B)
        _sync_and_announce("comm_external(emf_B) [after 2nd half B]")
        x.comm_external(runko.tools.comm_mode.emf_B)
        _sync_and_announce("comm_local(emf_B) [after 2nd half B]")
        x.comm_local(runko.tools.comm_mode.emf_B)

        # --- E push + edge BCs + add current ---
        _sync_and_announce("grid_push_e")
        x.grid_push_e()
        _sync_and_announce("grid_apply_edge_bcs(emf_E) [after E push]")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_E)
        _sync_and_announce("grid_add_current")
        x.grid_add_current()
        _sync_and_announce("grid_apply_edge_bcs(emf_E) [after add_current]")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_E)
        _sync_and_announce("comm_external(emf_E)")
        x.comm_external(runko.tools.comm_mode.emf_E)
        _sync_and_announce("comm_local(emf_E)")
        x.comm_local(runko.tools.comm_mode.emf_E)

        # --- advance wall position ---
        _sync_and_announce("prtcl_advance_reflector_walls")
        x.prtcl_advance_reflector_walls()

        # --- moving particle injector ---
        _sync_and_announce("injector.inject")
        injector.inject(simulation, [(0, pgen0), (1, pgen1)], ppc)

        # --- IO ---
        _sync_and_announce("io_average_kinetic_energy")
        x.io_average_kinetic_energy()
        _sync_and_announce("io_average_B_energy_density")
        x.io_average_B_energy_density()
        _sync_and_announce("io_average_E_energy_density")
        x.io_average_E_energy_density()
        _sync_and_announce("io_ram_usage")
        x.io_ram_usage()

        if simulation.lap % output_interval == 0:
            _sync_and_announce("io_emf_snapshot")
            x.io_emf_snapshot()
            #x.io_prtcl_snapshot()
            _sync_and_announce("io_spectra_snapshot")
            x.io_spectra_snapshot()
            _sync_and_announce("simulation.log_timer_statistics")
            simulation.log_timer_statistics()
            _sync_and_announce("simulation.reset_timers")
            simulation.reset_timers()

        _sync_and_announce("END OF LAP (returning from pic_simulation_step)")

    simulation.for_each_lap(pic_simulation_step)
    simulation.log_timer_statistics()
