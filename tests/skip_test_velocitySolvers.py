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




class MomentumInitialization(unittest.TestCase):

    def test_initialize(self):
        vsols1 = plasma.MomentumLagrangianSolver() 

    # here we just test that solver is able to take a step
    # NOTE: Correctness of the result is not tested here
    def test_ArbitrarySolve(self):
        self.Nx = 50
        self.Ny = 30
        self.Nz = 5

        self.params = Params()
        self.params.mins = [ -10.0, -10.0, -10.0 ]
        self.params.maxs = [  10.0,  10.0,  10.0 ]

        mesh = ptools.VeloMesh()
        mesh.Nblocks = [self.Nx, self.Ny, self.Nz]


        intp = ptools.BundleInterpolator4th()

        vsols = [ plasma.MomentumLagrangianSolver() ]
        for vsol in vsols:
            vsol.setInterpolator(intp)

            # TODO how to split cell tests and solver tests?
            #vsol.solve()





if __name__ == '__main__':
    unittest.main()
