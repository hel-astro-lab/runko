import unittest

import sys
sys.path.append('python')
import numpy as np

import plasmatools as ptools
import pyplasma as plasma


class Params:
    mins = None
    maxs = None
    lens = None



class Initialization(unittest.TestCase):

    i = 10
    j = 11
    o = 1

    def setUp(self):
        self.cell = plasma.VCell(self.i, self.j, self.o)

    # test that we can inherit from the corgi::Cell base class
    def test_inheritance(self):

        (i,j) = self.cell.index()
        self.assertEqual(i, 10)
        self.assertEqual(j, 11)

        self.cell.bark()
        



if __name__ == '__main__':
    unittest.main()
