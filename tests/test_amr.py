import unittest

import numpy as np
import sys

import pyplasma


class conf:

    outdir = "out"

    Nxv = 3
    Nyv = 4
    Nzv = 11

    xmin = -2.0
    ymin = -3.0
    zmin = -4.0

    xmax =  2.0
    ymax =  3.0
    zmax =  4.0


def gauss(ux,uy,uz):

    #return ux + uy + uz

    delgam = np.sqrt(1.0)
    mux = 0.0
    muy = 0.0
    muz = 0.0

    #f  = 1.0/np.sqrt(2.0*np.pi*delgam)
    f = 1.0
    f *= np.exp(-0.5*((ux - mux)**2)/delgam)
    f *= np.exp(-0.5*((uy - muy)**2)/delgam)
    f *= np.exp(-0.5*((uz - muz)**2)/delgam)

    return f



class Basics(unittest.TestCase):

    def setUp(self):

        self.m = pyplasma.AdaptiveMesh3D()
        self.m.resize( [conf.Nxv,  conf.Nyv,  conf.Nzv])
        self.m.set_min([conf.xmin, conf.ymin, conf.zmin])
        self.m.set_max([conf.xmax, conf.ymax, conf.zmax])


    def test_fill(self):

        for rfl in range(3):

            nx, ny, nz = self.m.get_size(rfl)
            comparison = np.zeros((nx,ny,nz))
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        x,y,z = self.m.get_center([i,j,k], rfl)
                        val = gauss(x,y,z)

                        self.m[i,j,k, rfl] =  val
                        comparison[i,j,k] = val


            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        val1 = self.m[i,j,k, rfl]
                        val2 = comparison[i,j,k]

                        self.assertAlmostEqual(val1, val2)


    def test_indexing(self):

        for rfl in range(3):
            nx, ny, nz = self.m.get_size(rfl)

            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):

                        #print(" ")

                        #set this in
                        indx1 = [i,j,k]
                        rfl1  = rfl
                        cid   = self.m.get_cell_from_indices(indx1, rfl1)
                        #print("indx: ", indx1, rfl1, " cid:", cid)


                        #and get this out
                        #print("test: ", indx1, " cid:", cid, " rfl1:", rfl1)
                        rfl2  = self.m.get_refinement_level(cid)
                        indx2 = self.m.get_indices(cid)
                        #print("indx: ", indx2, rfl2, " cid:", cid)

                        self.assertEqual(rfl1,  rfl2)
                        self.assertEqual(indx1, indx2)





if __name__ == '__main__':
    unittest.main()





