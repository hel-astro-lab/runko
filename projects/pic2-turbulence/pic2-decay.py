"""
Canonical pic simulation using runko with decaying fields.
"""


import pytools
import runko
import numpy as np
import itertools

if __name__ == "__main__":

    seed = 42
    np.random.seed(seed)

    conf = runko.Configuration(None)

    config.Nx = 4
    config.Ny = 4
    config.Nz = 4
    config.NxMesh = 20
    config.NyMesh = 20
    config.NzMesh = 20
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 0.45
    config.field_propagator = "FDTD2"
    config.q0 = -1
    config.m0 = 1
    config.q1 = 1
    config.m1 = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"


    # Problem specific configuration
    ppc = 2 # particles per cell (one particle type)
    oppc = 2 * ppc # overall particles per cell (all particle types)
    gamma = 1
    c_omp = 1
    omp = config.cfl / c

    config.q0 = -gamma * (omp**2.0) / (0.5 * oppc * (1.0 + config.m0 / config.m1))
    config.q1 *= self.q0

    config.m0 *= abs(config.q0)
    config.m1 *= abs(config.q1)

    delgam = 0.3 # temperature
    temp_ration = 1 # T_i / T_e
    sigma = 10 # temperature corrections
               # or magnetization number (omega_ce/omega_pe)^2, including gamma for inertia??

    delgam0 = delgam
    delgam1 = temp_ration * delgam0

    # No corrections; cold sigma
    binit_nc = sqrt(oppc * (config.cfl**2.0) * sigma * config.m0)
    # another approximation which is more accurate at \delta ~ 1
    gammath = 1.0 + (3.0/2.0) * delgam1

    binit_approx = sqrt(gammath * oppc * config.m0 * (config.cfl**2.0) * sigma)
    binit = binit_approx # NOTE: selecting this as our sigma definitions

    # Decaying setup:
    A0 = 0.8 * binit
    n_perp = 1
    n_par = 2
    ph1 = 2.0 * np.pi * np.random.rand(n_perp, n_perp, n_par)
    ph2 = 2.0 * np.pi * np.random.rand(n_perp, n_perp, n_par)
    ph3 = 2.0 * np.pi * np.random.rand(n_perp, n_perp, n_par)

    Lx = config.Nx * config.NxMesh / c_omp
    Ly = config.Ny * config.NyMesh / c_omp
    Lz = config.Nz * config.NzMesh / c_omp

    kx = 2.0 * np.pi / (c_omp * Lx)
    ky = 2.0 * np.pi / (c_omp * Ly)
    kz = 2.0 * np.pi / (c_omp * Lz)

    # Length of the wave mode vector
    def beta(n, m):
        return np.sqrt(8.0) / (np.sqrt(n**2 + m**2) * n_perp * np.sqrt(n_par))

    def B(x, y, z):
        Bx = 0
        By = 0
        Bz = binit

        I = range(1, n_perp + 1)
        for n, m in itertools.product(I, I):
            norm = beta(n ,m)
            for o in range(1, n_par):
                xmodx = np.sin(m * kx * x + ph1[n - 1, m - 1, o - 1])
                xmody = np.cos(n * ky * y + ph2[n - 1, m - 1, o - 1])

                ymodx = np.cos(m * kx * x + ph1[n - 1, m - 1, o - 1])
                ymody = np.sin(n * ky * y + ph2[n - 1, m - 1, o - 1])

                zmod = np.sin(o * kz * z + ph3[n - 1, m - 1, o - 1])

                Bx += norm * n * xmodx * xmody * zmod
                By -= norm * m * ymodx * ymody * zmod

        return A0 * Bx, A0 * By, Bz

    zero_field = lambda x, y, z: (0, 0, 0)

    def sample_vel(local_delgam):
        # velocity sampling from Maxwell-Juttner
        gamma = 0 # no bulk motion
        ux, uy, uz, _ = pytools.sample_boosted_maxwellian(local_delgam,
                                                          gamma,
                                                          direction=1,
                                                          dims=3)
        return ux, uy, uz

    def make_pgen(local_delgam):
        def f(x, y, z):
            loc = np.array((x, y, z))
            particles = []
            for _ in range(ppc):
                particles.append(runko.ParticleState(pos=loc + np.random.rand(),
                                                     vel=sample_vel(local_delgam)))
            return particles
        return f


    # TileGrid ctor:
    # - balances tiles based on conf (catepillar, hilbert)
    # - checks if restarts files are present for current config
    #   - if found initializes the tiles based on the restart files
    tile_grid = runko.TileGrid(conf)

    if not grid.initialized_from_restart_file():
        for idx in tile_grid.local_tile_indices():
            tile = runko.pypic.Tile(idx, conf) # Tile can deduce its coordinate extents from conf.
            tile.set_fields(zero_field, B, zero_field)
            tile.inject_to_each_cell(0, make_pgen(delgam0))
            tile.inject_to_each_cell(1, make_pgen(delgam1))
            tile_grid.add_tile(idx, tile)

    # Initializes the simulation and returns a handle to it:
    # - analyze and sync boundaries between mpi tasks/ranks/localities
    # - loads virtual tiles (is there benefit of doing this explicitly?)
    simulation = tile_grid.configure_simulation(conf)

    def sync_EB(_, comm, *_):
        EB = (runko.comm_mode.emf_E, runko.comm_mode.emf_B)
        comm.virtual_tile_sync(*EB)
        comm.pairwise_moore(*EB)

    simulation.prelude(sync_EB)

    def pic_simulation_step(tile, comm, *_):
        tile.push_half_b()
        comm.virtual_tile_sync(runko.comm_mode.emf_B)
        comm.pairwise_moore(runko.comm_mode.emf_B)

        tile.push_particles(0)
        tile.push_particles(1)
        comm.virtual_tile_sync(runko.comm_mode.pic_particle)
        comm.pairwise_moore(runko.comm_mode.pic_particle)

        tile.deposit_current()
        comm.virtual_tile_sync(runko.comm_mode.emf_J)
        comm.pairwise_moore(runko.comm_mode.emf_J_exchange)

        tile.push_half_b()
        comm.virtual_tile_sync(runko.comm_mode.emf_B)
        comm.pairwise_moore(runko.comm_mode.emf_B)

        tile.push_e()
        tile.add_J_to_E()
        comm.virtual_tile_sync(runko.comm_mode.emf_E)
        comm.pairwise_moore(runko.comm_mode.emf_E)

    simulation.for_each_lap(pic_simulation_step)
