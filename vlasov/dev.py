from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import sys, os

import pyplasmaDev as pdev


class conf:
    Nxv = 10
    Nyv = 10
    Nyz = 10

    xmin = -10.0
    ymin = -10.0
    zmin = -10.0

    xmax =  10.0
    ymax =  10.0
    zmax =  10.0



# set up the grid
m = pdev.AdaptiveMesh3D()
m.resize( [conf.Nxv,  conf.Nyv,  conf.Nzv])




    














