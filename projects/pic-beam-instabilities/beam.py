"""
Particle-in-cell simulation of relativistic beam instability.

Setup follows CfA Plasma Talks by Antoine Bret section 8 and 9
and its reference "Multidimensional electron beam-plasma instabilities
in the relativistic regime" by Bret, Gremillet and Dieckmann (2010).

Setup description: counter-streaming (macro-)electron beams `b` and `p`
over a bacground of fixed (macro-)protons. When protons are fixed
they do not generate any current and thus can be neglected.
Streams are aligend to x-direction.

Numerical values in code are prefixed with `num_` and
macro particles are denoted with `mp`.

Diluted two-stream (alpha << 1) is hard to produce with pic,
so atm this file only focuses on symmetric streaas (alpha = 1).

This files assumes that it is run with one MPI rank only.

Output is a csv file with following columns:

- normalized_time:
  t * omega_pe, where omega_pe is non-relativistic plasma frequency of beam p.

- B2_{par,perp}:
  average of B energy density {parallel,perpendicular} to the beam normalized
  by the relativistic plasma enthalphy density (gamma_b n_0 m_e c^2)

- E2_par, E2_perp:
  same as B2_{par,perp} but for E

- {two_stream,oblique,filamentation}_growth:
  represents the growth rate of the fastest growing mode for the three instabilities
  normalized such that initially the value is 1.
"""

import runko
import numpy as np


if __name__ == "__main__":
    logger = runko.runko_logger()

    config = runko.Configuration(None)

    # Input parameters:
    skin_depth_in_delta_x = 10
    num_n_0p = 32
    gamma_b = 3
    delgam_b = 1e-5 # delgam = k_b T / mc^2
    delgam_p = 1e-5
    config.Nt = 500
    config.cfl = 0.45
    field_output_interval = None # None means that the output is disabled.
    output_filename = "beam-instability" # .csv is automatically appended to the name
    enable_filter = False

    alpha = 1 # no other values supported.

    config.NxMesh = 320
    config.NyMesh = 80
    config.NzMesh = 6
    config.Nx = 1
    config.Ny = 1
    config.Nz = 1
    config.stride = 1
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.field_propagator_cfl_coeff = 1.02
    config.field_propagator = "FDTD2"
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st_atomic"
    config.current_filter = "binomial2"
    config.tile_partitioning = "hilbert_curve"


    num_mp_plasma_freq = config.cfl / skin_depth_in_delta_x

    num_n_0b = alpha * num_n_0p
    num_n = num_n_0p + num_n_0b

    v = lambda gamma: np.sqrt(1 - 1 / gamma**2)
    gamma = lambda v: 1.0 / np.sqrt(1 - v**2)

    v_b = v(gamma_b)
    v_p = -alpha * v_b
    gamma_p = gamma(v_p)

    gamma_mean = (alpha * gamma_b + gamma_p) / (1 + alpha)

    num_mp_charge = num_mp_plasma_freq**2 * gamma_mean / num_n

    config.q0 = -num_mp_charge
    config.m0 = 1

    # Classical (non-relativistic) plasma frequency for p.
    # Note that numerical mass and charge have the same value, so they cancel out.
    num_plasma_freq_p = np.sqrt(num_n_0p * abs(num_mp_charge))

    if alpha == 1:
        # delta_t from numerical growth rate and plasma freq cancel out.
        num_two_stream_growth_rate = num_plasma_freq_p / (2 * gamma_b**(3/2))
        num_oblique_growth_rate = num_plasma_freq_p * 0.5 * (3 / gamma_b)**(1/2)
        num_filamentation_growth_rate = num_plasma_freq_p * v_b * (2 / gamma_b)**(1/2)

    else:
        raise RuntimeError(f"This does not support alpha other than 1.")

    seed = 42
    logger.info(f"seed: {seed}")
    logger.info(f"alpha: {alpha}")
    growth_rate_log = f"{num_two_stream_growth_rate:.03E} "
    growth_rate_log += f"{num_oblique_growth_rate:.03E} "
    growth_rate_log += f"{num_filamentation_growth_rate:.03E}"
    logger.info(f"growth rates (two-stream, oblique, filamentation): {growth_rate_log}")
    logger.info(f"num_n_0b: {num_n_0b}")
    logger.info(f"num_n_0p: {num_n_0p}")
    logger.info(f"gamma (b, p, mean): {gamma_b:.03E}, {gamma_p:.03E}, {gamma_mean:.03E}")
    logger.info(f"q0, m0: {config.q0:.03E}, {config.m0:.03E}")
    logger.info(f"delgam (b, p): {delgam_b:.03E}, {delgam_p:.03E}")

    zero_field = lambda x, y, z: np.zeros_like(x)

    def pgen_beam(num_n, gamma, delgam, dir: str):
        def pgen(x, y, z):
            # Use the same seed such that the beam ions are generated on top of each other.
            rng = np.random.default_rng(seed=seed)

            x, y, z = np.repeat(x, num_n), np.repeat(y, num_n), np.repeat(z, num_n)

            N = len(x)
            dx, dy, dz = rng.random(N), rng.random(N), rng.random(N)
            pos = x + dx, y + dy, z + dz

            vel = runko.sample_boosted_juttner_synge(N, delgam, Gamma=gamma, direction=dir)
            return runko.pic.threeD.ParticleStateBatch(pos=pos, vel=vel)

        return pgen

    tile_grid = runko.TileGrid(config)

    if not tile_grid.initialized_from_restart_file():
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            tile.batch_set_EBJ(zero_field, zero_field, zero_field,
                               zero_field, zero_field, zero_field,
                               zero_field, zero_field, zero_field)

            logger.info("injecting b...")
            tile.batch_inject_to_cells(0, pgen_beam(num_n_0b, gamma_b, delgam_b, "x"))
            logger.info("injecting p...")
            tile.batch_inject_to_cells(0, pgen_beam(num_n_0p, gamma_p, delgam_p, "-x"))

            tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(config)

    laps = []
    B2_parallel_avgs = []
    E2_parallel_avgs = []
    B2_perpendicular_avgs = []
    E2_perpendicular_avgs = []

    def store_data():
        B2_parallel_sum = 0
        E2_parallel_sum = 0
        B2_perpendicular_sum = 0
        E2_perpendicular_sum = 0

        N = 0
        for tile in simulation.local_tiles():
            (Ex, Ey, Ez), (Bx, By, Bz), _ = tile.get_EBJ()

            N += Ex.size

            E2_parallel_sum += np.sum(Ex**2)
            E2_perpendicular_sum += np.sum(Ey**2 + Ez**2)
            B2_parallel_sum += np.sum(Bx**2)
            B2_perpendicular_sum += np.sum(By**2 + Bz**2)

        laps.append(simulation.lap)
        E2_parallel_avgs.append(E2_parallel_sum / N)
        E2_perpendicular_avgs.append(E2_perpendicular_sum / N)
        B2_parallel_avgs.append(B2_parallel_sum / N)
        B2_perpendicular_avgs.append(B2_perpendicular_sum / N)

    def sync_EB(tile, comm, io):
        EB = (runko.tools.comm_mode.emf_E, runko.tools.comm_mode.emf_B)
        comm.virtual_tile_sync(*EB)
        comm.pairwise_moore(*EB)

    simulation.prelude(sync_EB)

    def pic_simulation_step(tile, comm, io):
        store_data()

        tile.push_half_b()
        comm.virtual_tile_sync(runko.tools.comm_mode.emf_B)
        comm.pairwise_moore(runko.tools.comm_mode.emf_B)

        tile.push_particles()
        comm.virtual_tile_sync(runko.tools.comm_mode.pic_particle)
        comm.pairwise_moore(runko.tools.comm_mode.pic_particle)

        if simulation.lap % 5 == 0:
            tile.sort_particles()

        tile.deposit_current()
        comm.virtual_tile_sync(runko.tools.comm_mode.emf_J)
        comm.pairwise_moore(runko.tools.comm_mode.emf_J_exchange)
        comm.virtual_tile_sync(runko.tools.comm_mode.emf_J)
        comm.pairwise_moore(runko.tools.comm_mode.emf_J)
        if enable_filter:
            tile.filter_current()
            comm.virtual_tile_sync(runko.tools.comm_mode.emf_J)
            comm.pairwise_moore(runko.tools.comm_mode.emf_J)
            tile.filter_current()
            tile.filter_current()


        tile.push_half_b()
        comm.virtual_tile_sync(runko.tools.comm_mode.emf_B)
        comm.pairwise_moore(runko.tools.comm_mode.emf_B)

        tile.push_e()
        tile.subtract_J_from_E()
        comm.virtual_tile_sync(runko.tools.comm_mode.emf_E)
        comm.pairwise_moore(runko.tools.comm_mode.emf_E)

        if field_output_interval and simulation.lap % field_output_interval == 0:
            io.emf_snapshot()

        if simulation.lap % 20 == 0:
            simulation.log_timer_statistics()


    simulation.for_each_lap(pic_simulation_step)
    store_data()
    simulation.log_timer_statistics()

    laps = np.array(laps)
    normalized_t = laps * num_plasma_freq_p

    enthalphy_coeff = 2 * gamma_b * config.cfl**2 * abs(num_mp_charge) * num_n

    B2_parallel = np.array(B2_parallel_avgs) / enthalphy_coeff
    B2_perpendicular = np.array(B2_perpendicular_avgs) / enthalphy_coeff

    E2_parallel = np.array(E2_parallel_avgs) / enthalphy_coeff
    E2_perpendicular = np.array(E2_perpendicular_avgs) / enthalphy_coeff

    two_stream_growth = np.exp(num_two_stream_growth_rate * laps)
    oblique_growth = np.exp(num_oblique_growth_rate * laps)
    filamentation_growth = np.exp(num_filamentation_growth_rate * laps)

    table = np.column_stack((normalized_t,
                             B2_parallel,
                             B2_perpendicular,
                             E2_parallel,
                             E2_perpendicular,
                             two_stream_growth,
                             oblique_growth,
                             filamentation_growth))

    header = "normalized_time, "
    header += "B2_par, "
    header += "B2_perp, "
    header += "E2_par, "
    header += "E2_perp, "
    header += "two_stream_growth, "
    header += "oblique_growth, "
    header += "Filamentation_growth"

    np.savetxt(output_filename + ".csv",
               table,
               delimiter=', ',
               header=header)
