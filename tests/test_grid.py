import unittest
import numpy as np

import pyplasmabox




class Initialization(unittest.TestCase):


    def test_node1D(self):

        node = pyplasmabox.vlv.Grid1D(3)
        node.setGridLims(0.0, 1.0)

        #c = pyplasmabox.fields.Tile1D(10)
        c = pyplasmabox.vlv.vlvTile1D(10)
        node.addTile(c, (1,) ) 

        c2 = node.getTile(1)

        (i) = c2.index
        self.assertEqual( 1, i )


    def test_node2D(self):

        node = pyplasmabox.vlv.Grid2D(3,3)
        node.setGridLims(0.0, 1.0, 10.0, 20.0)

        # note: we also try different tile family here
        c = pyplasmabox.fields.Tile2D(10, 10)
        node.addTile(c, (1,2) ) 

        c2 = node.getTile(1,2)

        (i,j) = c2.index
        self.assertEqual( 1, i )
        self.assertEqual( 2, j )


