import numpy as np
from numpy import sqrt, pi

import pytools

# Field initialization
def insert_em_fields(grid, conf):
    from numpy import arctan2, sin, cos, sqrt, pi

    Lx = conf.Nx*conf.NxMesh
    Ly = conf.Ny*conf.NyMesh
    Lz = conf.Nz*conf.NzMesh

    b0 = conf.binit   # equilibrium magnetic field
    zeta = conf.zeta  # perturbation amplitude

    # maximum perpendicular and parallel wavenumbers
    kperp = 1.0*pi/Lx #half of the sine only
    kpar  = 2.0*pi/Lz

    #phase shift
    om0 = 0.0 

    #amplitude; 
    A = b0*zeta/2.0

    # position of the centers as Stagger objects
    pkg_loc1 = conf.pkg_loc1
    pkg_loc2 = conf.pkg_loc2 

    #-------------------------------------------------- 
    #ell = Lx/2. #maximum size
    ell = conf.ell

    for cid in grid.get_tile_ids():
        tile = grid.get_tile(cid)
        gs = tile.get_grids(0)

        if conf.twoD:
            ii, jj = tile.index
            kk = 0
        elif conf.threeD:
            ii, jj, kk = tile.index

        # insert values into Yee lattices
        for n in range(conf.NzMesh):
            for m in range(conf.NyMesh):
                for l in range(conf.NxMesh):

                    # get global coordinates
                    iglob, jglob, kglob = pytools.ind2loc((ii, jj, kk), (l, m, n), conf)
                    loc = pytools.Stagger(iglob, jglob, kglob)

                    # distance from the center of the packets
                    d1 = pkg_loc1 - loc
                    d2 = pkg_loc2 - loc

                    x = iglob
                    y = jglob
                    z = kglob 

                    # add stationary guide field
                    gs.bz[l, m, n] = b0

                    bpkg1 = {}
                    bpkg2 = {}

                    # build exact initial amplitude for different staggered grid locations
                    for st in [ "rh", "bx", "by", "bz", "ex", "ey", "ez",]:

                        # spherical distance
                        r1 = np.sqrt( d1.at(st).x ** 2 + d1.at(st).y ** 2 + d1.at(st).z ** 2 )
                        r2 = np.sqrt( d2.at(st).x ** 2 + d2.at(st).y ** 2 + d2.at(st).z ** 2 )

                        envelope1 = np.cos(pi*r1/ell/2.0)
                        envelope2 = np.cos(pi*r2/ell/2.0)

                        # cut after r > ell
                        if r1 > 1.10*ell or envelope1 < 0.0:
                            envelope1 = 0.0
                        if r2 > 1.10*ell or envelope2 < 0.0:
                            envelope2 = 0.0

                        # volume normalized amplitdue
                        #bpkg1[st] = b0*(pi * ell) * envelope
                        #bpkg2[st] = b0*(pi * ell) * envelope

                        # peak normalized amplitdue
                        bpkg1[st] = 2.0*b0*zeta*envelope1/ell
                        bpkg2[st] = 2.0*b0*zeta*envelope2/ell


                    #first pkg
                    gs.bx[l, m, n] += -bpkg1["bx"] * d1.at("bx").y
                    gs.by[l, m, n] += +bpkg1["by"] * d1.at("by").x

                    wdir = +1.0*conf.beta #direction of wave
                    gs.ex[l, m, n] += wdir*bpkg1["ex"] * d1.at("ex").x
                    gs.ey[l, m, n] += wdir*bpkg1["ey"] * d1.at("ey").y

                    #second pkg
                    if conf.two_wave:
                        if conf.two_wave_reversed:
                            rdir = -1.0 #direction of twist; antialigned
                        else:
                            rdir = +1.0 #direction of twist; aligned

                        gs.bx[l, m, n] += -rdir*bpkg2["bx"] * d2.at("bx").y
                        gs.by[l, m, n] += +rdir*bpkg2["by"] * d2.at("by").x

                        wdir = -1.0*conf.beta #direction of wave
                        gs.ex[l, m, n] += wdir*rdir*bpkg2["ex"] * d2.at("ex").x
                        gs.ey[l, m, n] += wdir*rdir*bpkg2["ey"] * d2.at("ey").y

    return



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
        delgam = conf.delgam_e

    # positrons/ions/second species
    elif ispcs == 1:
        delgam = conf.delgam_i

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand(1)
    yy = xloc[1] + np.random.rand(1)
    zz = xloc[2] + np.random.rand(1)

    # velocity sampling
    gamma = 0.0
    direction = +1
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



