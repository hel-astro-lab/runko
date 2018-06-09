import unittest

import sys
sys.path.append('python')
import numpy as np

import corgi
import pyplasma as plasma
import pypic 

sys.path.append('pic')
from pic import loadCells
from pic import inject
from pic import spatialLoc


from visualize import saveVisz

try:
    import matplotlib.pyplot as plt
except:
    pass



def const_field(x, y, z):
    return 1.0


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(node, conf, ffunc):

    Lx  = conf.Nx*conf.NxMesh #XXX scaled length
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            c = node.getCellPtr(i,j)
            yee = c.getYee(0)

            for l in range(conf.NxMesh):
                for m in range(conf.NyMesh):
                    for n in range(conf.NzMesh):

                        #get x_i+1/2 (Yee lattice so rho_i)
                        xloc0 = spatialLoc(node, (i,j), (l,  m,n), conf)
                        xloc1 = spatialLoc(node, (i,j), (l+1,m,n), conf)

                        xmid = 0.5*(xloc0[0] + xloc1[0])
                        ymid = 0.5*(xloc0[1] + xloc1[1])
                        zmid = 0.5*(xloc0[2] + xloc1[2])

                        val = ffunc(xmid, ymid, zmid)

                        yee.ex[l,m,n] = val
                        yee.ey[l,m,n] = val+1.0
                        yee.ez[l,m,n] = val+2.0

                        yee.bx[l,m,n] = val+3.0
                        yee.by[l,m,n] = val+4.0
                        yee.bz[l,m,n] = val+5.0

                        yee.jx[l,m,n] = val
                        yee.jy[l,m,n] = val
                        yee.jz[l,m,n] = val


# basic Conf file/class for PiC simulation testing
class Conf:

    Nx = 1
    Ny = 1
    Nz = 1

    NxMesh = 5
    NyMesh = 5
    NzMesh = 1

    xmin = 0.0
    xmax = 10.0

    ymin = 0.0
    ymax = 10.0

    cfl = 0.45
    c_omp = 10.0
    ppc = 1

    dx = 1.0
    dy = 1.0
    dz = 1.0

    me = 1
    mi = 1

    Nspecies = 1

    outdir = "out"

    #def __init__(self):
    #    print("initialized...")

    #update bounding box sizes
    def update_bbox(self):
        self.xmin = 0.0
        self.xmax = self.Nx*self.NxMesh

        self.ymin = 0.0
        self.ymax = self.Ny*self.NyMesh


# 2d -> 1d
def flatten(arr):
    return arr.flatten().tolist()

# 1d -> 2d
def reshape(vec, nx, ny):
    return np.reshape(vec, (nx, ny))

class Filters(unittest.TestCase):



    def test_init(self):

        NxMesh = 3
        NyMesh = 3

        #internal filter size is 3x Nx/y/zMesh
        NxF = NxMesh*3
        NyF = NyMesh*3

        flt = pypic.Filter(NxMesh, NyMesh)

        kernel = np.random.rand(NxF, NyF)
        image  = np.random.rand(NxF, NyF)

        flt.set_kernel( flatten(kernel) )
        flt.set_image(  flatten(image) )

        k2   = flt.get_kernel()
        img2 = flt.get_image()

        k2   = reshape(k2,   NxF, NyF)
        img2 = reshape(img2, NxF, NyF)

        for i in range(NxF):
            for j in range(NyF):
                self.assertEqual(kernel[i,j], k2[i,j])
                self.assertEqual(image[i,j], img2[i,j])
                

    def skip_test_kernel_init(self):

        NxMesh = 3
        NyMesh = 3

        #internal filter size is 3x Nx/y/zMesh
        NxF = NxMesh*3
        NyF = NyMesh*3

        flt = pypic.Filter(NxMesh, NyMesh)

        flt.init_kernel()
        #kernel = np.random.rand(NxF, NyF)
        #flt.set_kernel( flatten(kernel) )

        #image  = np.random.rand(NxF, NyF)
        #flt.set_image(  flatten(image) )

        k2   = reshape( flt.get_kernel(), NxF, NyF)
        #img2 = reshape( flt.get_image() , NxF, NyF)
                
        print()
        print(k2)

        flt.fft_kernel()
        k3   = reshape( flt.get_kernel(), NxF, NyF)

        print()
        print(k3)



    # FFT transform image forward and then backward to see if we get the same result back
    # NOTE: floating-point conversion if not exact so we need some error tolerance
    def test_fft_backandforth(self):

        NxMesh = 3
        NyMesh = 3

        #internal filter size is 3x Nx/y/zMesh
        NxF = NxMesh*3
        NyF = NyMesh*3

        flt = pypic.Filter(NxMesh, NyMesh)

        #flt.init_kernel()
        #kernel = np.random.rand(NxF, NyF)
        #flt.set_kernel( flatten(kernel) )

        #create and set image
        image  = np.random.rand(NxF, NyF)
        flt.set_image(  flatten(image) )
        img1 = reshape( flt.get_image() , NxF, NyF)
        #print("orig")
        #print(img1)

        flt.fft_image_forward()
        img2 = reshape( flt.get_image() , NxF, NyF)
        #print("fft forward")
        #print(img2)

        flt.fft_image_backward()
        img3 = reshape( flt.get_image() , NxF, NyF)
        #print("fft backward")
        #print(img3)

        for i in range(NxF):
            for j in range(NyF):
                self.assertAlmostEqual(img1[i,j], img3[i,j], places=6) 




    def test_filters_in_action(self):

        conf = Conf()
        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 1
        conf.NxMesh = 5
        conf.NyMesh = 5
        conf.NzMesh = 1

        node = plasma.Grid(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadCells(node, conf)
        insert_em(node, conf, const_field)
        #inject(node, filler_no_velocity, conf) #injecting plasma particles




