from mpi4py import MPI
import unittest
import numpy as np

import pycorgi
import pyplasmabox




class Initialization(unittest.TestCase):


    def test_node1D(self):

        node = pycorgi.oneD.Node(3,1,1)
        node.setGridLims(0.0, 1.0)

        #c = pyplasmabox.fields.Tile1D(10)
        c = pyplasmabox.vlv.oneD.Tile(10,1,1)
        node.addTile(c, (1,) ) 

        c2 = node.getTile(1)

        (i,) = c2.index
        self.assertEqual( 1, i )


    def test_node2D(self):

        node = pycorgi.twoD.Node(3,3,1)
        node.setGridLims(0.0, 1.0, 10.0, 20.0)

        # note: we also try different tile family here
        c = pyplasmabox.fields.twoD.Tile(10, 10, 1)
        node.addTile(c, (1,2) ) 

        c2 = node.getTile(1,2)

        (i,j) = c2.index
        self.assertEqual( 1, i )
        self.assertEqual( 2, j )


