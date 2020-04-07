import numpy as np
import pytools  # runko python tools


# Prtcl velocity (and location modulation inside cell)
#
# NOTE: Cell extents are from xloc to xloc + 1, i.e., dx = 1 in all dims.
#       Therefore, typically we use position as x0 + RUnif[0,1).
#
# NOTE: Default injector changes odd ispcs's loc to (ispcs-1)'s prtcl loc.
#       This means that positrons/ions are on top of electrons to guarantee
#       charge conservation (because zero charge density initially).
#
def velocity_profile(xloc, ispcs, conf):

    # electrons
    if ispcs == 0:
        delgam = conf.delgam  # * np.abs(conf.mi / conf.me) * conf.temp_ratio
        direction = -1

    # positrons/ions/second species
    elif ispcs == 1:
        delgam = conf.delgam
        direction = +1

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    zz = xloc[2] + np.random.rand(1)

    # velocity sampling
    gamma = conf.gamma

    ux, uy, uz, uu = pytools.sample_boosted_maxwellian(
        delgam, gamma, direction=direction, dims=3
    )

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0


# number of prtcls of species 'ispcs' to be added to cell at location 'xloc'
#
# NOTE: Plasma frequency is adjusted based on conf.ppc (and prtcl charge conf.qe/qi
#       and mass conf.me/mi are computed accordingly) so modulate density in respect
#       to conf.ppc only.
#
def density_profile(xloc, ispcs, conf):

    # uniform plasma with default n_0 number density
    return conf.ppc


# Field initialization
def insert_em_fields(grid, conf):

    # into radians
    btheta = conf.btheta / 180.0 * np.pi
    bphi = conf.bphi / 180.0 * np.pi
    beta = conf.beta

    for tile in pytools.tiles_all(grid):
        yee = tile.get_yee(0)

        ii,jj,kk = tile.index if conf.threeD else (*tile.index, 0)

        # insert values into Yee lattices; includes halos from -3 to n+3
        for n in range(-3, conf.NzMesh + 3):
            for m in range(-3, conf.NyMesh + 3):
                for l in range(-3, conf.NxMesh + 3):
                    # get global coordinates
                    iglob, jglob, kglob = pytools.ind2loc((ii, jj, kk), (l, m, n), conf)
                    r = np.sqrt(iglob ** 2 + jglob ** 2 + kglob ** 2)

                    yee.bx[l, m, n] = conf.binit * np.cos(bphi)
                    yee.by[l, m, n] = conf.binit * np.sin(bphi) * np.sin(btheta)
                    yee.bz[l, m, n] = conf.binit * np.sin(bphi) * np.cos(btheta)

                    yee.ex[l, m, n] = 0.0
                    yee.ey[l, m, n] = -beta * yee.bz[l, m, n]
                    yee.ez[l, m, n] = beta * yee.by[l, m, n]
    return

