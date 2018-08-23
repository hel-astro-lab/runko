import unittest
import numpy as np

import pyplasmabox as plasma


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
        self.tile = plasma.VlasovTile()

    # test that we can inherit from the corgi::Tile base class
    def test_inheritance(self):

        (i,j) = self.tile.index
        self.assertEqual(i, 10)
        self.assertEqual(j, 11)

        



if __name__ == '__main__':
    unittest.main()
