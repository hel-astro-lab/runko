from mpi4py import MPI
import unittest

import os

import numpy as np
import pycorgi
import pyrunko
import pytools
import h5py

#from visualize import get_yee
#from visualize import getYee2D
#from combine_files import combine_tiles
#import injector
#from read_mesh import TileInfo
#from read_mesh import get_mesh
#import initialize as init



# combine tiles inside grid together into one array
def combine_tiles(ff, fvar, conf, isp=None ):

    arr = np.zeros((conf.Nx*conf.NxMesh, 
                    conf.Ny*conf.NyMesh, 
                    conf.Nz*conf.NzMesh))

    f = h5py.File(ff,'r')
    for dset in f:

        if not(isp==None):
            if not(f[dset]['ispcs'][()] == isp):
                continue

        i = f[dset]['i'][()]
        j = f[dset]['j'][()]
        k = f[dset]['k'][()]

        NxMesh = f[dset]['Nx'][()]
        NyMesh = f[dset]['Ny'][()]
        NzMesh = f[dset]['Nz'][()]
    
        ii = int( i*NxMesh )
        jj = int( j*NyMesh )
        kk = int( k*NzMesh )

        tile = f[dset][fvar][()]
        tile = np.reshape(tile, (NzMesh, NyMesh, NxMesh))

        for s in range(NzMesh):
            for r in range(NyMesh):
                for q in range(NxMesh):
                    arr[ii+q, jj+r, kk+s] = tile[s,r,q]

    return arr


class Conf:

    Nx = 3
    Ny = 3
    Nz = 1

    oneD = False
    twoD = False
    threeD = False

    NxMesh = 5
    NyMesh = 5
    NzMesh = 5

    xmin = 0.0
    xmax = 1.0

    ymin = 0.0
    ymax = 1.0

    zmin = 0.0
    zmax = 1.0

    me = -1.0
    mi =  1.0

    qe = 1.0
    qi = 1.0

    #def __init__(self):
    #    print("initialized...")


def density_profile(xloc, ispcs, conf):
    return conf.ppc

#load tiles into each grid
def loadTiles1D(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #if n.get_mpi_grid(i) == n.rank:
            c = pyrunko.fields.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            n.add_tile(c, (i,) ) 


#load tiles into each grid
def loadTiles2D(n, conf):
    for i in range(n.get_Nx()):
        for j in range(n.get_Ny()):
            #if n.get_mpi_grid(i,j) == n.rank:
            c = pyrunko.fields.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            n.add_tile(c, (i,j) ) 


# create similar reference array
def fill_ref(grid, conf):
    data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 9))

    # lets put floats into Yee lattice
    val = 1.0
    Nx = grid.get_Nx()
    Ny = grid.get_Ny()
    Nz = grid.get_Nz()

    NxM = conf.NxMesh
    NyM = conf.NyMesh
    NzM = conf.NzMesh

    #print("filling ref")
    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            for k in range(grid.get_Nz()):
                #if n.get_mpi_grid(i,j) == n.rank:
                if True:
                    for q in range(conf.NxMesh):
                        for r in range(conf.NyMesh):
                            for s in range(conf.NzMesh):
                                #print(i,j,k,q,r,s, "val=",val)
                                # insert 1/val to get more complex floats (instead of ints only)
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 0] = 1.0/val + 1.0 
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 1] = 1.0/val + 2.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 2] = 1.0/val + 3.0

                                data[i*NxM + q, j*NyM + r, k*NzM + s, 3] = 1.0/val + 4.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 4] = 1.0/val + 5.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 5] = 1.0/val + 6.0

                                data[i*NxM + q, j*NyM + r, k*NzM + s, 6] = 1.0/val + 7.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 7] = 1.0/val + 8.0
                                data[i*NxM + q, j*NyM + r, k*NzM + s, 8] = 1.0/val + 9.0
                                val += 1
    return data


# fill Yee mesh with values
def fill_yee(grid, data, conf):

    Nx = grid.get_Nx()
    Ny = grid.get_Ny()
    Nz = grid.get_Nz()

    NxM = conf.NxMesh
    NyM = conf.NyMesh
    NzM = conf.NzMesh

    # lets put ref array into Yee lattice
    #print("filling yee")
    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):
            for k in range(grid.get_Nz()):
                #if n.get_mpi_grid(i,j) == n.rank:
                if True:
                    c = grid.get_tile(i,j,k)
                    yee = c.get_yee(0)

                    for q in range(conf.NxMesh):
                        for r in range(conf.NyMesh):
                            for s in range(conf.NzMesh):
                                yee.ex[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 0] 
                                yee.ey[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 1]
                                yee.ez[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 2]

                                yee.bx[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 3]
                                yee.by[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 4]
                                yee.bz[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 5]

                                yee.jx[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 6]
                                yee.jy[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 7]
                                yee.jz[q,r,s] = data[i*NxM + q, j*NyM + r, k*NzM + s, 8]


# Generic function to fill the velocity mesh
def filler(xloc, uloc, ispcs, conf):
    x = xloc[0]
    y = xloc[1]
    z = xloc[2] 

    ux = uloc[0]
    uy = uloc[1]
    uz = uloc[2] 

    return x + y + z + ux + uy + uz + ispcs

    #physical version
    mux = 0.0
    muy = 0.0
    muz = 0.0
    delgam = conf.delgam

    #electrons
    if ispcs == 0:
        delgam  = conf.delgam * np.abs(conf.mi / conf.me) * conf.temperature_ratio
        mux = conf.ub_e
        muy = 0.0
        muz = 0.0
    #positrons/ions/second species
    if ispcs == 1:
        delgam  = conf.delgam
        mux = conf.ub_i
        muy = 0.0
        muz = 0.0

    #plasma reaction
    omp = conf.cfl*conf.dx
    n0 = (omp**2.0)/conf.Nspecies

    #velocity perturbation
    Lx = conf.Nx*conf.NxMesh*conf.dx
    kmode = conf.modes
    mux_noise = conf.beta*np.cos(2.0*np.pi*kmode*x/Lx) * (Lx/(2.0*np.pi*kmode))

    #Classical Maxwellian
    f  = n0*(1.0/(2.0*np.pi*delgam))**(0.5)
    f *= np.exp(-0.5*((ux - mux - mux_noise)**2.0)/(delgam))
    return f


class IO(unittest.TestCase):

    def test_write_fields1D(self):

        ##################################################
        # write
        conf = Conf()
        conf.oneD = True

        conf.Nx = 3
        conf.Ny = 1
        conf.Nz = 1
        conf.NxMesh = 2
        conf.NyMesh = 1
        conf.NzMesh = 1 
        conf.outdir = "io_test_1D/"

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        grid = pycorgi.oneD.Grid(conf.Nx, conf.Ny)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles1D(grid, conf)

        ref = fill_ref(grid, conf)
        fill_yee(grid, ref, conf)

        pyrunko.fields.oneD.write_yee(grid, 0, conf.outdir)

        ##################################################
        # read using analysis tools

        arrs = combine_tiles(conf.outdir+"fields-0_0.h5", "ex", conf)

        Nx = grid.get_Nx()
        Ny = grid.get_Ny()
        Nz = grid.get_Nz()

        NxM = conf.NxMesh
        NyM = conf.NyMesh
        NzM = conf.NzMesh
        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                for k in range(grid.get_Nz()):
                    #if n.get_mpi_grid(i,j) == n.rank:
                    if True:
                        for q in range(conf.NxMesh):
                            for r in range(conf.NyMesh):
                                for s in range(conf.NzMesh):
                                    self.assertAlmostEqual( arrs[i*NxM + q, j*NyM + r, k*NzM + s],  
                                                             ref[i*NxM + q, j*NyM + r, k*NzM + s, 0], places=6)

        ##################################################
        # test reading back
        node2 = pycorgi.oneD.Grid(conf.Nx, conf.Ny)
        node2.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles1D(node2, conf)

        pyrunko.fields.oneD.read_yee(node2, 0, "io_test_1D")

        yee1 = pytools.visualize.get_yee(grid,  conf)
        yee2 = pytools.visualize.get_yee(node2, conf)

        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                for k in range(grid.get_Nz()):
                    #if n.get_mpi_grid(i,j) == n.rank:
                    if True:
                        for q in range(conf.NxMesh):
                            for r in range(conf.NyMesh):
                                for s in range(conf.NzMesh):
                                    for m in ["jx", "jy", "jz", "ex", "ey", "ez", "bx", "by", "bz", "rho"]:

                                        self.assertAlmostEqual( 
                                                yee1[m][i*NxM + q],
                                                yee2[m][i*NxM + q],
                                                places=6)



    def test_write_fields2D(self):

        ##################################################
        # write

        conf = Conf()
        conf.twoD = True

        conf.Nx = 3
        conf.Ny = 4
        conf.Nz = 1
        conf.NxMesh = 5
        conf.NyMesh = 6
        conf.NzMesh = 1 
        conf.outdir = "io_test_2D/"

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles2D(grid, conf)

        ref = fill_ref(grid, conf)
        fill_yee(grid, ref, conf)

        pyrunko.fields.twoD.write_yee(grid, 0, conf.outdir)
        
        ##################################################
        # read using analysis tools
        arrs = combine_tiles(conf.outdir+"fields-0_0.h5", "ex", conf)

        Nx = grid.get_Nx()
        Ny = grid.get_Ny()
        Nz = grid.get_Nz()

        NxM = conf.NxMesh
        NyM = conf.NyMesh
        NzM = conf.NzMesh
        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                for k in range(grid.get_Nz()):
                    #if n.get_mpi_grid(i,j) == n.rank:
                    if True:
                        for q in range(conf.NxMesh):
                            for r in range(conf.NyMesh):
                                for s in range(conf.NzMesh):
                                    self.assertAlmostEqual( arrs[i*NxM + q, j*NyM + r, k*NzM + s],  
                                                             ref[i*NxM + q, j*NyM + r, k*NzM + s, 0], places=6)


        ##################################################
        # test reading back
        node2 = pycorgi.twoD.Grid(conf.Nx, conf.Ny)
        node2.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles2D(node2, conf)

        pyrunko.fields.twoD.read_yee(node2, 0, "io_test_2D")

        yee1 = pytools.visualize.get_yee_2D(grid,  conf)
        yee2 = pytools.visualize.get_yee_2D(node2, conf)

        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                for k in range(grid.get_Nz()):
                    #if n.get_mpi_grid(i,j) == n.rank:
                    if True:
                        for q in range(conf.NxMesh):
                            for r in range(conf.NyMesh):
                                for s in range(conf.NzMesh):
                                    for m in ["jx", "jy", "jz", "ex", "ey", "ez", "bx", "by", "bz", "rho"]:

                                        self.assertAlmostEqual( 
                                                yee1[m][i*NxM + q, j*NyM + r],
                                                yee2[m][i*NxM + q, j*NyM + r],
                                                places=6)


    # compare two AdaptiveMesh3D objects and assert their equality
    def compareMeshes(self, vm, ref):
        cells = vm.get_cells(True)
        refcells = ref.get_cells(True)

        #metainfo
        self.assertEqual( vm.length, ref.length )
        self.assertEqual( vm.maximum_refinement_level, ref.maximum_refinement_level )
        self.assertEqual( vm.top_refinement_level, ref.top_refinement_level )
        self.assertEqual( len(cells), len(refcells) )

        for cid in cells:
            #refinement level
            rfl1 = vm.get_refinement_level(cid)
            rfl2 = ref.get_refinement_level(cid)
            self.assertEqual(rfl1, rfl2)

            #indices
            [ii1,jj1,kk1] = vm.get_indices(cid)
            [ii2,jj2,kk2] =ref.get_indices(cid)
            self.assertEqual(ii1, ii2)
            self.assertEqual(jj1, jj2)
            self.assertEqual(kk1, kk2)
        
            #value
            self.assertEqual( vm[ii1,jj1,kk1,rfl1], ref[ii2,jj2,kk2,rfl2] )

            #center
            xx1,yy1,zz1 = vm.get_center([ii1,jj1,kk1], rfl1)
            xx2,yy2,zz2 =ref.get_center([ii2,jj2,kk2], rfl2)
            self.assertEqual(xx1, xx2)
            self.assertEqual(yy1, yy2)
            self.assertEqual(zz1, zz2)



    # Complicated restart from file with heterogeneous (boundary) tiles
    def skip_test_restart(self):

        conf = Conf()
        conf.twoD = True

        conf.Nx = 5
        conf.Ny = 1
        conf.Nz = 1
        conf.NxMesh = 5
        conf.NyMesh = 1
        conf.NzMesh = 1 
        conf.outdir = "io_test_restart/"

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        grid = pycorgi.oneD.Grid(conf.Nx, conf.Ny)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        pytools.vlv.loadTiles(grid, conf)
    
        #load boundaries
        for i in [0, grid.get_Nx()-1]:
            for j in range(grid.get_Ny()):
                if i == 0:
                    c = plasma.Tile_outflow_L(conf.NxMesh, conf.NyMesh, conf.NzMesh)

                    #set left (top) wall location
                    c.fld1 = 1
                    c.fld2 = 1

                else:
                    c = plasma.Tile_outflow_R(conf.NxMesh, conf.NyMesh, conf.NzMesh)

                    #set right (bottom) wall location
                    l = conf.NxMesh-2
                    iglob, jglob, kglob = globalIndx( (i,j), (l,0,0), conf)

                    c.fld1 = iglob
                    c.fld2 = iglob
                
                pytools.vlv.initialize_tile(c, (i,j), grid, conf)

                #add it to the grid
                grid.add_tile(c, (i,)) 



    def test_write_pic2D(self):

        def test_filler(xloc, ispcs, conf):
        
            xx = xloc[0] 
            yy = xloc[1] 

            #electrons
            if ispcs == 0:
                zz = 0.1
        
            #positrons/ions/second species
            if ispcs == 1:
                zz = 0.2
        
            ux = xx*100.0
            uy = yy*1000.0
            uz = xx*yy*10000.0
        
            x0 = [xx, yy, zz]
            u0 = [ux, uy, uz]

            return x0, u0



        ##################################################
        # write

        conf = Conf()
        conf.twoD = True

        conf.Nx = 3
        conf.Ny = 4
        conf.Nz = 1
        conf.NxMesh = 5
        conf.NyMesh = 6
        conf.NzMesh = 1 
        conf.outdir = "io_test_2D/"
        conf.ppc = 1
        conf.Nspecies = 2
        conf.Nspecies_test = 0

        #tmp non-needed variables
        conf.omp = 1
        conf.gamma_e = 0.0
        conf.me = 1
        conf.mi = 1
        conf.cfl = 1.0
        conf.c_omp = 1.0


        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        for i in range(grid.get_Nx()):
            for j in range(grid.get_Ny()):
                for k in range(grid.get_Nz()):
                    c = pyrunko.pic.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                    pytools.pic.initialize_tile(c, (i, j, k), grid, conf)
                    grid.add_tile(c, (i,j)) 

        pytools.pic.inject(grid, test_filler, density_profile, conf)

        pyrunko.pic.twoD.write_particles(grid, 0, conf.outdir)

        # TODO: read with h5py



        # TODO: read with internal read tool
        node2 = pycorgi.twoD.Grid(conf.Nx, conf.Ny)
        node2.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        for i in range(node2.get_Nx()):
            for j in range(node2.get_Ny()):
                for k in range(node2.get_Nz()):
                    c = pyrunko.pic.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                    pytools.pic.initialize_tile(c, (i, j,k), node2, conf)
                    node2.add_tile(c, (i,j)) 

        pyrunko.pic.twoD.read_particles(node2, 0, conf.outdir)

        #assert
        for i in range(node2.get_Nx()):
            for j in range(node2.get_Ny()):
                if node2.get_mpi_grid(i,j) == node2.rank():
                    cid    = node2.id(i,j)
                    c      = node2.get_tile(cid) #get cell ptr

                    c_ref  = grid.get_tile(cid) #get cell ptr

                    for ispcs in range(conf.Nspecies):
                        container1 = c.get_container(ispcs)
                        container2 = c_ref.get_container(ispcs)

                        #TODO: assert content

                        xxs1 = container1.loc(0)
                        yys1 = container1.loc(1)
                        zzs1 = container1.loc(2)
                        vxs1 = container1.vel(0)
                        vys1 = container1.vel(1)
                        vzs1 = container1.vel(2)
                        wgs1 = container1.wgt()

                        xxs2 = container2.loc(0)
                        yys2 = container2.loc(1)
                        zzs2 = container2.loc(2)
                        vxs2 = container2.vel(0)
                        vys2 = container2.vel(1)
                        vzs2 = container2.vel(2)
                        wgs2 = container2.wgt()

                        nprtcls = conf.NxMesh*conf.NyMesh*conf.NzMesh*conf.ppc

                        self.assertEqual(len(xxs1), len(xxs2))
                        self.assertEqual(len(yys1), len(yys2))
                        self.assertEqual(len(zzs1), len(zzs2))
                        
                        self.assertEqual(len(vxs1), len(vxs2))
                        self.assertEqual(len(vys1), len(vys2))
                        self.assertEqual(len(vzs1), len(vzs2))

                        self.assertEqual(len(wgs1), len(wgs2))

                        for n in range(nprtcls):
                            self.assertAlmostEqual(xxs1[n], xxs2[n], places=6)
                            self.assertAlmostEqual(yys1[n], yys2[n], places=6)
                            self.assertAlmostEqual(zzs1[n], zzs2[n], places=6)

                            self.assertAlmostEqual(vxs1[n], vxs2[n], places=6)
                            self.assertAlmostEqual(vys1[n], vys2[n], places=6)
                            self.assertAlmostEqual(vzs1[n], vzs2[n], places=6)

                            self.assertAlmostEqual(wgs1[n], wgs2[n], places=6)












