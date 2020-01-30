# -*- coding: utf-8 -*- 

import numpy as np

import pycorgi
import pyrunko



def balance_mpi(n, comm_size=None):

    nx = n.get_Nx()
    ny = n.get_Ny()
    nz = n.get_Nz()

    if nz > 1:
        D = 3
    elif ny > 1:
        D = 2
    elif nx > 1:
        D = 1

    print("Dimensionality of grid: {}".format(D))

    if D == 2:
        return balance_mpi_2D(n, comm_size=comm_size)
    elif D == 3:
        return balance_mpi_3D(n, comm_size=comm_size)



# load nodes using 2D Hilbert curve
def balance_mpi_2D(n, comm_size=None):

    if n.master:  # only master initializes; then sends
        if comm_size == None:
            comm_size = n.size()

        nx = n.get_Nx()
        ny = n.get_Ny()

        m0 = np.log2(nx)
        m1 = np.log2(ny)

        if not (m0.is_integer()):
            raise ValueError("Nx is not power of 2 (i.e. 2^m)")

        if not (m1.is_integer()):
            raise ValueError("Ny is not power of 2 (i.e. 2^m)")

        # print('Generating hilbert with 2^{} {}'.format(m0,m1))
        hgen = pyrunko.tools.twoD.HilbertGen(np.int(m0), np.int(m1))

        igrid = np.zeros((nx, ny), np.int64)
        grid = np.zeros((nx, ny))  # , np.int64)

        for i in range(nx):
            for j in range(ny):
                grid[i, j] = hgen.hindex(i, j)

        # print(grid)
        hmin, hmax = np.min(grid), np.max(grid)
        for i in range(nx):
            for j in range(ny):
                igrid[i, j] = np.floor(comm_size * grid[i, j, k] / (hmax + 1))

        # check that nodes get about same work load
        # y = np.bincount(igrid.flatten())
        # ii = np.nonzero(y)[0]
        # print(list(zip(ii,y[ii])))

        # print("grid:")
        for i in range(nx):
            for j in range(ny):
                val = igrid[i, j]
                n.set_mpi_grid(i, j, val)
                    # print("({},{}) = {}".format(i,j,val))

    # broadcast calculation to everybody
    n.bcast_mpi_grid()


# load nodes using 3D Hilbert curve
def balance_mpi_3D(n, comm_size=None):

    if n.master:  # only master initializes; then sends
        if comm_size == None:
            comm_size = n.size()

        nx = n.get_Nx()
        ny = n.get_Ny()
        nz = n.get_Nz()

        m0 = np.log2(nx)
        m1 = np.log2(ny)
        m2 = np.log2(nz)

        if not (m0.is_integer()):
            raise ValueError("Nx is not power of 2 (i.e. 2^m)")

        if not (m1.is_integer()):
            raise ValueError("Ny is not power of 2 (i.e. 2^m)")

        if not (m2.is_integer()):
            raise ValueError("Nz is not power of 2 (i.e. 2^m)")

        # print('Generating hilbert with 2^{} {}'.format(m0,m1))
        hgen = pyrunko.tools.threeD.HilbertGen(np.int(m0), np.int(m1), np.int(m2))

        igrid = np.zeros((nx, ny, nz), np.int64)
        grid = np.zeros((nx, ny, nz))  # , np.int64)

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    grid[i, j, k] = hgen.hindex(i, j, k)

        # print(grid)
        hmin, hmax = np.min(grid), np.max(grid)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    igrid[i, j, k] = np.floor(comm_size * grid[i, j, k] / (hmax + 1))

        # check that nodes get about same work load
        # y = np.bincount(igrid.flatten())
        # ii = np.nonzero(y)[0]
        # print(list(zip(ii,y[ii])))

        # print("grid:")
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    val = igrid[i, j, k]
                    n.set_mpi_grid(i, j, k, val)
                    # print("({},{}) = {}".format(i,j,val))

    # broadcast calculation to everybody
    n.bcast_mpi_grid()
