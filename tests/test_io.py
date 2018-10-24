from mpi4py import MPI
import unittest

import os

import numpy as np
import pycorgi
import pyplasmabox
import h5py

from visualize import getYee
from visualize import getYee2D
from combine_files import combine_tiles

import injector
from read_mesh import TileInfo
from read_mesh import get_mesh

import initialize as init

class Conf:

    Nx = 3
    Ny = 3
    Nz = 1

    NxMesh = 5
    NyMesh = 5
    NzMesh = 5

    xmin = 0.0
    xmax = 1.0

    ymin = 0.0
    ymax = 1.0

    me = -1.0
    mi =  1.0

    #def __init__(self):
    #    print("initialized...")



#load tiles into each node
def loadTiles1D(n, conf):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #if n.getMpiGrid(i) == n.rank:
            c = pyplasmabox.fields.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            n.addTile(c, (i,) ) 


#load tiles into each node
def loadTiles2D(n, conf):
    for i in range(n.getNx()):
        for j in range(n.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            c = pyplasmabox.fields.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
            n.addTile(c, (i,j) ) 


# create similar reference array
def fill_ref(node, conf):
    data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 9))

    # lets put floats into Yee lattice
    val = 1.0
    Nx = node.getNx()
    Ny = node.getNy()
    Nz = node.getNz()

    NxM = conf.NxMesh
    NyM = conf.NyMesh
    NzM = conf.NzMesh

    #print("filling ref")
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            for k in range(node.getNz()):
                #if n.getMpiGrid(i,j) == n.rank:
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
def fill_yee(node, data, conf):

    Nx = node.getNx()
    Ny = node.getNy()
    Nz = node.getNz()

    NxM = conf.NxMesh
    NyM = conf.NyMesh
    NzM = conf.NzMesh

    # lets put ref array into Yee lattice
    #print("filling yee")
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            for k in range(node.getNz()):
                #if n.getMpiGrid(i,j) == n.rank:
                if True:
                    c = node.getTile(i,j,k)
                    yee = c.getYee(0)

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
        print("working in fields1D: ", os.getcwd())

        conf = Conf()
        conf.Nx = 3
        conf.Ny = 1
        conf.Nz = 1
        conf.NxMesh = 2
        conf.NyMesh = 1
        conf.NzMesh = 1 
        conf.outdir = "io_test_1D/"

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        node = pycorgi.oneD.Node(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles1D(node, conf)

        ref = fill_ref(node, conf)
        fill_yee(node, ref, conf)

        pyplasmabox.vlv.oneD.writeYee(node, 0, conf.outdir)
        
        print("wrote hdf5")
        ##################################################
        # read using analysis tools

        arrs = combine_tiles(conf.outdir+"fields-0_0.h5", "ex", conf)

        Nx = node.getNx()
        Ny = node.getNy()
        Nz = node.getNz()

        NxM = conf.NxMesh
        NyM = conf.NyMesh
        NzM = conf.NzMesh
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                for k in range(node.getNz()):
                    #if n.getMpiGrid(i,j) == n.rank:
                    if True:
                        for q in range(conf.NxMesh):
                            for r in range(conf.NyMesh):
                                for s in range(conf.NzMesh):
                                    self.assertAlmostEqual( arrs[i*NxM + q, j*NyM + r, k*NzM + s],  
                                                             ref[i*NxM + q, j*NyM + r, k*NzM + s, 0], places=6)

        print("read using analysis tools")
        ##################################################
        # test reading back
        node2 = pycorgi.oneD.Node(conf.Nx, conf.Ny)
        node2.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles1D(node2, conf)

        print("readinfg using ezh5")
        pyplasmabox.vlv.oneD.readYee(node2, 0, "io_test_1D")
        print("read using c++ tools")

        yee1 = getYee(node,  conf)
        yee2 = getYee(node2, conf)

        for i in range(node.getNx()):
            for j in range(node.getNy()):
                for k in range(node.getNz()):
                    #if n.getMpiGrid(i,j) == n.rank:
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
        print("working in fields2D: ", os.getcwd())

        conf = Conf()
        conf.Nx = 3
        conf.Ny = 4
        conf.Nz = 1
        conf.NxMesh = 5
        conf.NyMesh = 6
        conf.NzMesh = 1 
        conf.outdir = "io_test_2D/"

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        node = pycorgi.twoD.Node(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        loadTiles2D(node, conf)

        ref = fill_ref(node, conf)
        fill_yee(node, ref, conf)

        pyplasmabox.vlv.twoD.writeYee(node, 0, conf.outdir)
        
        ##################################################
        # read using analysis tools
        arrs = combine_tiles(conf.outdir+"fields-0_0.h5", "ex", conf)

        Nx = node.getNx()
        Ny = node.getNy()
        Nz = node.getNz()

        NxM = conf.NxMesh
        NyM = conf.NyMesh
        NzM = conf.NzMesh
        for i in range(node.getNx()):
            for j in range(node.getNy()):
                for k in range(node.getNz()):
                    #if n.getMpiGrid(i,j) == n.rank:
                    if True:
                        for q in range(conf.NxMesh):
                            for r in range(conf.NyMesh):
                                for s in range(conf.NzMesh):
                                    self.assertAlmostEqual( arrs[i*NxM + q, j*NyM + r, k*NzM + s],  
                                                             ref[i*NxM + q, j*NyM + r, k*NzM + s, 0], places=6)


        ##################################################
        # test reading back
        node2 = pycorgi.twoD.Node(conf.Nx, conf.Ny)
        node2.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles2D(node2, conf)

        pyplasmabox.vlv.twoD.readYee(node2, 0, "io_test_2D")

        yee1 = getYee2D(node,  conf)
        yee2 = getYee2D(node2, conf)

        for i in range(node.getNx()):
            for j in range(node.getNy()):
                for k in range(node.getNz()):
                    #if n.getMpiGrid(i,j) == n.rank:
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


    def test_write_Mesh3V(self):
        print("working in mesh3V: ", os.getcwd())

        ##################################################
        # write
        conf = Conf()
        conf.Nx = 2
        conf.Ny = 1
        conf.Nz = 1
        conf.NxMesh = 3
        conf.NyMesh = 1
        conf.NzMesh = 1 
        conf.Nspecies = 2
        conf.outdir = "io_test_mesh/"

        conf.dx  = 1.0
        conf.dy  = 1.0
        conf.dz  = 1.0
        conf.Nvx = 3
        conf.Nvy = 5
        conf.Nvz = 8

        conf.vxmin =  -1.0
        conf.vymin =  -2.0
        conf.vzmin =  -3.0
        conf.vxmax =   4.0
        conf.vymax =   5.0
        conf.vzmax =   6.0
        conf.refinement_level= 0
        conf.clip= True
        conf.clipThreshold= 1.0e-5

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        node = pycorgi.oneD.Node(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        for i in range(node.getNx()):
            for j in range(node.getNy()):
                #if n.getMpiGrid(i) == n.rank:
                c = pyplasmabox.vlv.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                node.addTile(c, (i,) ) 
        injector.inject(node, filler, conf ) #injecting plasma

        pyplasmabox.vlv.oneD.writeMesh(node, 0, conf.outdir)

        ##################################################
        # read using analysis tools
        fname = conf.outdir + "meshes-0_0.h5"
        f = h5py.File(fname,'r')

        for i in range(node.getNx()):
            for j in range(node.getNy()):
                for k in range(node.getNz()):
                    c = node.getTile(i,j,k)

                    #if n.getMpiGrid(i,j) == n.rank:
                    if True:
                        for ispcs in range(conf.Nspecies):
                            block = c.getPlasmaSpecies(0, ispcs)

                            for q in range(conf.NxMesh):
                                for r in range(conf.NyMesh):
                                    for s in range(conf.NzMesh):
                                        tinfo = TileInfo()
                                        tinfo.i = i
                                        tinfo.j = j
                                        tinfo.k = k
                                        tinfo.q = q
                                        tinfo.r = r
                                        tinfo.s = s
                                        tinfo.ispcs = ispcs

                                        #now assert 
                                        vm = get_mesh(f, tinfo)
                                        ref = block[q,r,s]
                                        self.compareMeshes(vm, ref)



        ##################################################
        # read back

        node2 = pycorgi.oneD.Node(conf.Nx, conf.Ny)
        node2.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        for i in range(node2.getNx()):
            for j in range(node2.getNy()):
                #if n.getMpiGrid(i) == n.rank:
                c = pyplasmabox.vlv.oneD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                node2.addTile(c, (i,) ) 
        injector.inject(node2, injector.empty_filler, conf, empty=True) #injecting empty meshes

        #pyplasmabox.vlv.oneD.writeMesh(node2, 1, conf.outdir)
        pyplasmabox.vlv.oneD.readMesh(node2,  0, "io_test_mesh")
        #pyplasmabox.vlv.oneD.writeMesh(node2, 2, conf.outdir)

        for i in range(node2.getNx()):
            for j in range(node2.getNy()):
                for k in range(node2.getNz()):
                    c1 = node.getTile(i,j,k)
                    c2 = node2.getTile(i,j,k)

                    #if n.getMpiGrid(i,j) == n.rank:
                    if True:
                        for ispcs in range(conf.Nspecies):
                            block1 = c1.getPlasmaSpecies(0, ispcs)
                            block2 = c2.getPlasmaSpecies(0, ispcs)

                            for q in range(conf.NxMesh):
                                for r in range(conf.NyMesh):
                                    for s in range(conf.NzMesh):
                                        tinfo = TileInfo()
                                        tinfo.i = i
                                        tinfo.j = j
                                        tinfo.k = k
                                        tinfo.q = q
                                        tinfo.r = r
                                        tinfo.s = s
                                        tinfo.ispcs = ispcs

                                        #now assert 
                                        vm1 = block1[q,r,s]
                                        vm2 = block2[q,r,s]
                                        self.compareMeshes(vm1, vm2)



    # Complicated restart from file with heterogeneous (boundary) tiles
    def skip_test_restart(self):

        conf = Conf()
        conf.Nx = 5
        conf.Ny = 1
        conf.Nz = 1
        conf.NxMesh = 5
        conf.NyMesh = 1
        conf.NzMesh = 1 
        conf.outdir = "io_test_restart/"

        if not os.path.exists( conf.outdir ):
            os.makedirs(conf.outdir)

        node = pycorgi.oneD.Node(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        init.loadTiles(node, conf)
    
        #load boundaries
        for i in [0, node.getNx()-1]:
            for j in range(node.getNy()):
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
                
                init.initialize_tile(c, i,j, node, conf)

                #add it to the node
                node.addTile(c, (i,)) 









