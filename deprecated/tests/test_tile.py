from mpi4py import MPI
import unittest
import numpy as np

import pyrunko


class Params:
    mins = None
    maxs = None
    lens = None



class Initialization(unittest.TestCase):

    i = 10
    j = 11
    o = 1
    Nx = 10
    Ny = 5

    NxMesh = 3
    NyMesh = 3
    NzMesh = 1

    def setUp(self):
        self.tile = pyrunko.emf.twoD.Tile(self.NxMesh, self.NyMesh, self.NzMesh)
        self.tile.index = (self.i, self.j)

    # test that we can inherit from the corgi::Tile base class
    def test_inheritance(self):

        (i,j) = self.tile.index
        self.assertEqual(i, self.i)
        self.assertEqual(j, self.j)

        



if __name__ == '__main__':
    unittest.main()
