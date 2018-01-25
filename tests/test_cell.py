import unittest

import sys
sys.path.append('python')
import numpy as np

import pyplasmaDev
import pyplasma as plasma


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

    def setUp(self):
        self.cell = plasma.VlasovCell(self.i, self.j, self.o, 
                                      self.Nx, self.Ny, 
                                      self.NxMesh, self.NyMesh)

    # test that we can inherit from the corgi::Cell base class
    def test_inheritance(self):

        (i,j) = self.cell.index()
        self.assertEqual(i, 10)
        self.assertEqual(j, 11)

        



if __name__ == '__main__':
    unittest.main()
