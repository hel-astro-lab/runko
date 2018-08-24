import unittest
import numpy as np

import pyplasmabox


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
        self.tile = pyplasmabox.fields.Tile2D(self.NxMesh, self.NyMesh)
        self.tile.index = (self.i, self.j)

    # test that we can inherit from the corgi::Tile base class
    def test_inheritance(self):

        (i,j) = self.tile.index
        self.assertEqual(i, self.i)
        self.assertEqual(j, self.j)

        



if __name__ == '__main__':
    unittest.main()
