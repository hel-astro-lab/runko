from mpi4py import MPI
import unittest
import numpy as np

import pycorgi
import pyrunko




class Initialization(unittest.TestCase):

    #TODO add 1d pic
    #def test_node1D(self):
    #    grid = pycorgi.oneD.Grid(3,1,1)
    #    grid.set_grid_lims(0.0, 1.0)

    #    #c = pyrunko.emf.Tile1D(10)
    #    c = pyrunko.pic.oneD.Tile(10,1,1)
    #    grid.add_tile(c, (1,) ) 

    #    c2 = grid.get_tile(1)

    #    (i,) = c2.index
    #    self.assertEqual( 1, i )


    def test_node2D(self):

        grid = pycorgi.twoD.Grid(3,3,1)
        grid.set_grid_lims(0.0, 1.0, 10.0, 20.0)

        # note: we also try different tile family here
        c = pyrunko.emf.twoD.Tile(10, 10, 1)
        grid.add_tile(c, (1,2) ) 

        c2 = grid.get_tile(1,2)

        (i,j) = c2.index
        self.assertEqual( 1, i )
        self.assertEqual( 2, j )


