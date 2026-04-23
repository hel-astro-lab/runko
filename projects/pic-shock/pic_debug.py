"""
Particle-in-cell simulation of a collisionless shock — DEBUG copy.

Sibling of pic.py with full per-kernel profiling written to
<outdir>/ram-usage/per-step-0.csv (rank 0 only). Per lap, for every
kernel / comm / edge_bc / filter call in the step, records:

  - VmRSS delta (kB)            — host memory change
  - wall time (us)              — time.perf_counter_ns delta
  - Python traced memory at lap start/end (bytes, tracemalloc)
  - GPU device memory at lap start/end (kB, None->-1 on CPU backend)

The filter loop unrolls into a probe per pass + per comm pair, so
n_filter_passes changes the column count; the header is written
after the first lap. Delete this file (or discard the
v5-sublap-ram-trace branch) when the culprit is fixed.
"""

import runko
import numpy as np
import argparse
import logging
import os
import csv
import atexit
import time
import tracemalloc
from runko.ram_usage import get_rss_kB, get_gpu_mem_kB


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

    # --------------------------------------------------
    # Sub-lap profiling trace (DEBUG) — rank 0 only.
    #
    # Header is deferred to the first lap's write (so filter-loop
    # probes adapt to n_filter_passes automatically). Row layout:
    #   lap, rss_start_kB, py_mem_start_B, gpu_mem_start_kB,
    #   d_rss_<tag>...   (kB deltas, one per kernel call),
    #   t_us_<tag>...    (microsecond deltas, same order),
    #   rss_end_kB, py_mem_end_B, gpu_mem_end_kB
    _trace_f = None
    _trace_writer = None
    _trace_header_written = False
    _trace_tag_order: list[str] = []
    if runko.on_main_rank():
        from runko.auto_outdir import resolve_outdir as _resolve_outdir
        _trace_dir = os.path.join(_resolve_outdir(conf), "ram-usage")
        os.makedirs(_trace_dir, exist_ok=True)
        _trace_f = open(os.path.join(_trace_dir, "per-step-0.csv"), "w",
                        buffering=1)
        _trace_writer = csv.writer(_trace_f)
        atexit.register(_trace_f.close)
        tracemalloc.start()
        logger.info(f"sub-lap profile -> {_trace_dir}/per-step-0.csv "
                    "(rss/time per kernel + py/gpu at lap boundaries)")

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
    # Main PIC loop

    def pic_simulation_step(x):
        # --- DEBUG: per-kernel RSS + wall-time probes ---
        on_main = runko.on_main_rank()
        marks: list[tuple[str, int, int]] = []  # (tag, rss_kB, t_ns)
        def mark(tag):
            if on_main:
                marks.append((tag, get_rss_kB(), time.perf_counter_ns()))

        mark("start")
        if on_main:
            _py_mem_start, _ = tracemalloc.get_traced_memory()
            _gpu_kb = get_gpu_mem_kB()
            _gpu_mem_start = -1 if _gpu_kb is None else _gpu_kb

        # --- half B push #1 + B edge BCs + comm ---
        x.grid_push_half_b()                                       ; mark("push_half_b_1")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_B)         ; mark("edgebc_B_1")
        x.comm_external(runko.tools.comm_mode.emf_B)               ; mark("extB_1")
        x.comm_local(runko.tools.comm_mode.emf_B)                  ; mark("locB_1")

        # --- particle push + reflect + pack + exchange ---
        x.prtcl_push()                                             ; mark("prtcl_push")
        x.prtcl_reflect_particles()                                ; mark("prtcl_reflect")
        x.prtcl_pack_outgoing()                                    ; mark("prtcl_pack")
        x.comm_external(runko.tools.comm_mode.pic_particle)        ; mark("extP")
        x.comm_local(runko.tools.comm_mode.pic_particle)           ; mark("locP")

        if simulation.lap % 5 == 0:
            x.prtcl_sort()
        mark("prtcl_sort")

        # --- current deposit + J exchange + J edge_bc ---
        x.prtcl_deposit_current()                                  ; mark("deposit")
        x.comm_external(runko.tools.comm_mode.emf_J)               ; mark("extJ_1")
        x.comm_local(runko.tools.comm_mode.emf_J_exchange)         ; mark("locJ_exch")
        x.comm_external(runko.tools.comm_mode.emf_J)               ; mark("extJ_2")
        x.comm_local(runko.tools.comm_mode.emf_J)                  ; mark("locJ")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_J)         ; mark("edgebc_J_1")

        # --- current filter loop (unrolled into per-pass probes) ---
        # Optimal: communicate every halo_size (=3) passes. With
        # n_filter_passes > 3, two per-loop comms fire at i=3, 6, ...
        for i in range(n_filter_passes):
            if i > 0 and i % 3 == 0:
                x.comm_external(runko.tools.comm_mode.emf_J)       ; mark(f"filt_ext_{i}")
                x.comm_local(runko.tools.comm_mode.emf_J)          ; mark(f"filt_loc_{i}")
            x.grid_filter_current()                                ; mark(f"filt_{i}")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_J)         ; mark("edgebc_J_2")

        # --- half B push #2 + B edge BCs + comm ---
        x.grid_push_half_b()                                       ; mark("push_half_b_2")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_B)         ; mark("edgebc_B_2")
        x.comm_external(runko.tools.comm_mode.emf_B)               ; mark("extB_2")
        x.comm_local(runko.tools.comm_mode.emf_B)                  ; mark("locB_2")

        # --- E push + edge_bc + add_current + comm ---
        x.grid_push_e()                                            ; mark("push_e")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_E)         ; mark("edgebc_E_1")
        x.grid_add_current()                                       ; mark("add_current")
        x.grid_apply_edge_bcs(runko.tools.comm_mode.emf_E)         ; mark("edgebc_E_2")
        x.comm_external(runko.tools.comm_mode.emf_E)               ; mark("extE")
        x.comm_local(runko.tools.comm_mode.emf_E)                  ; mark("locE")

        # --- advance walls + inject ---
        x.prtcl_advance_reflector_walls()                          ; mark("advance_walls")
        injector.inject(simulation, [(0, pgen0), (1, pgen1)], ppc) ; mark("inject")

        # --- IO per lap ---
        x.io_average_kinetic_energy()                              ; mark("io_avg_ke")
        x.io_average_B_energy_density()                            ; mark("io_avg_B")
        x.io_average_E_energy_density()                            ; mark("io_avg_E")
        x.io_ram_usage()                                           ; mark("io_ram")

        # --- snapshot (every output_interval) ---
        if simulation.lap % output_interval == 0:
            x.io_emf_snapshot()
            #x.io_prtcl_snapshot()
            x.io_spectra_snapshot()
            simulation.log_timer_statistics()
            simulation.reset_timers()
        mark("io_snapshot")

        # --- DEBUG: flush one row ---
        if on_main and _trace_writer is not None:
            global _trace_header_written, _trace_tag_order
            tags = [m[0] for m in marks]
            if not _trace_header_written:
                _trace_tag_order = tags
                header = ["lap", "rss_start_kB",
                          "py_mem_start_B", "gpu_mem_start_kB"]
                header += [f"d_rss_{t}" for t in tags[1:]]
                header += [f"t_us_{t}" for t in tags[1:]]
                header += ["rss_end_kB", "py_mem_end_B", "gpu_mem_end_kB"]
                _trace_writer.writerow(header)
                _trace_header_written = True
            elif tags != _trace_tag_order:
                # If mark order drifts (e.g. n_filter_passes changes
                # mid-run), fall back to a sparse write — shouldn't
                # happen in practice.
                pass

            py_mem_end, _ = tracemalloc.get_traced_memory()
            _gpu_kb = get_gpu_mem_kB()
            gpu_mem_end = -1 if _gpu_kb is None else _gpu_kb

            rss_start = marks[0][1]
            rss_end = marks[-1][1]
            rss_deltas = [marks[i][1] - marks[i - 1][1]
                          for i in range(1, len(marks))]
            t_deltas_us = [(marks[i][2] - marks[i - 1][2]) // 1000
                           for i in range(1, len(marks))]
            row = [simulation.lap, rss_start,
                   _py_mem_start, _gpu_mem_start]
            row += rss_deltas
            row += t_deltas_us
            row += [rss_end, py_mem_end, gpu_mem_end]
            _trace_writer.writerow(row)

    simulation.for_each_lap(pic_simulation_step)
    simulation.log_timer_statistics()
