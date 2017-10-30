import unittest

import sys
sys.path.append('python')
import numpy as np

import plasmatools as plasma


class Params:
    mins = None
    maxs = None
    lens = None


def cellID2index(cellID, dvs):

    cellID -= 1

    k = np.int(  ( cellID / (dvs[0] * dvs[1]) ) )
    j = np.int(  ( cellID / dvs[0] ) % dvs[1] )
    i = np.int(  cellID % dvs[0] )

    return (i,j,k)


def populate_mesh( mesh ):

    for k in range(mesh.Nblocks[2]):
        for j in range(mesh.Nblocks[1]):
            for i in range(mesh.Nblocks[0]):
                cid = mesh.getBlockID([i,j,k])
                (x,y,z) = mesh.getCenter( cid )

                fval = physical_vel(x,y,z)
                mesh[i,j,k] = [fval, fval, fval, fval]
                #print "({},{},{}) = {}".format(i,j,k,fval)


#physical "real" distribution to compare against
def physical_vel(x,y,z):

    mux = 1.0
    muy = 2.0
    muz = 3.0
    sigmax = 2.0
    sigmay = 3.0
    sigmaz = 4.0

    vx = np.exp(-(x-mux)**2 / sigmax**2 )
    vy = np.exp(-(y-muy)**2 / sigmay**2 )
    #vz = np.exp(-(z-muz)**2 / sigmaz**2 )
    vz = 1.0

    return vx*vy*vz
    #return 0.5



class Basics(unittest.TestCase):

    def setUp(self):

        self.Nx = 50
        self.Ny = 30
        self.Nz = 5

        self.params = Params()
        self.params.mins = [ -10.0, -10.0, -10.0 ]
        self.params.maxs = [  10.0,  10.0,  10.0 ]

        self.mesh = plasma.VeloMesh()
        self.mesh.Nblocks = [self.Nx, self.Ny, self.Nz]


    def test_zFill(self):
        self.mesh.zFill( self.params.mins, self.params.maxs )

        self.assertEqual( self.mesh.number_of_blocks, self.Nx*self.Ny*self.Nz )

    def test_indices(self):

        tests = [ [1,1,0], [2,2,0], [3,3,0] ]
        for test in tests:

            #get cell id and compare to internal python formula
            cid = self.mesh.getBlockID( test )
            indx = cellID2index( cid , self.mesh.Nblocks )
            ref_indx = self.mesh.getIndices( cid )

            #check i, j, k indices
            self.assertEqual(indx[0], ref_indx[0])
            self.assertEqual(indx[1], ref_indx[1])
            self.assertEqual(indx[2], ref_indx[2])






class Data(unittest.TestCase):

    def setUp(self):

        self.Nx = 50
        self.Ny = 30
        self.Nz = 5

        self.params = Params()
        self.params.mins = [ -10.0, -10.0, -10.0 ]
        self.params.maxs = [  10.0,  10.0,  10.0 ]

        self.mesh = plasma.VeloMesh()
        self.mesh.Nblocks = [self.Nx, self.Ny, self.Nz]
        self.mesh.zFill( self.params.mins, self.params.maxs )

    def test_meshData(self):
        data = [1.0, 2.0, 3.0, 4.5]

        #test cellid indexing
        self.mesh[1] = data
        self.assertEqual( self.mesh[1], data )

        #test ijk indexing
        self.mesh[2,2,0] = data
        self.assertEqual( self.mesh[2,2,0], data )






if __name__ == '__main__':
    unittest.main()
