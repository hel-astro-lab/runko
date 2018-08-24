import sys
sys.path.append('python')
import numpy as np
from math import floor, ceil
from scipy.signal import convolve2d
from scipy.signal import convolve

import pycorgi
import pyplasmabox.pic as pypic

from initialize_pic import loadTiles
from initialize_pic import spatialLoc
from injector_pic import inject

from visualize import saveVisz

try:
    import matplotlib.pyplot as plt
except:
    pass



def const_field(x, y, z):
    return 1.0

def linear_ramp(x,y,z):
    return x + y + z


# insert initial electromagnetic setup (or solve Poisson eq)
def insert_em(node, conf, ffunc):

    Lx  = conf.Nx*conf.NxMesh #XXX scaled length
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            c = node.getTilePtr(i,j)
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

                        yee.jx[l,m,n] = val
                        yee.jy[l,m,n] = val+1.0
                        yee.jz[l,m,n] = val+2.0


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

    gamma_e = 0.0
    gamma_i = 0.0

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
    #return np.transpose(arr).flatten().tolist()

# 1d -> 2d
def reshape(vec, nx, ny):
    return np.reshape(vec, (nx, ny))
    #return np.transpose( np.reshape(vec, (nx, ny)) )
    #return np.reshape(np.flipud(vec), (nx, ny)) 

def flip(arr):
    return np.fliplr(np.flipud(arr))

#def fftshift1d(i, w):
#    if i >= 0:
#        return i
#    else:
#        return i - w

def fftshift1d(i, w):
    #return i + np.int(np.ceil(1.0*w/2.0))
    return i + np.int(np.floor(1.0*w/2.0))-1
    #return i + w/2



# fft shift 2D array
def fftshift(farr):
    nx, ny = np.shape(farr)
    arr = np.empty_like(farr)
    arr[:] = 0.2
    xc1 = np.int( np.floor(1.0*nx/2.0))
    yc1 = np.int( np.floor(1.0*ny/2.0))
    xc2 = np.int( np.ceil(1.0*nx/2.0) )
    yc2 = np.int( np.ceil(1.0*ny/2.0) )
    print(xc1, yc1)
    print(xc2, yc2)

    print("center (0,0) = ({},{})".format(fftshift1d(0,nx), fftshift1d(0,ny)))

    for i in np.arange(-xc1+1, xc2+1, 1):
        for j in np.arange(-yc1+1, yc2+1, 1):
            iff = fftshift1d(i, nx)
            jff = fftshift1d(j, ny)
            #print("ij=({},{}) i2j2=({},{})".format(i,j,iff,jff))
            if (iff < 0 or jff < 0):
                print("ERROR")

            #print("i={} j={} i2={} j2={}".format(i,j,iff,jff))
            arr[iff,jff] = farr[i,j]
    return arr


def sinc2d(kernel, X, Y):
    nx, ny = np.shape(kernel)
    xc = nx/3
    yc = ny/3

    for i in range(-xc, xc, 1):
        for j in range(-yc, yc, 1):
            iff = fftshift1d(i, nx)
            jff = fftshift1d(j, ny)
            arr[iff,jff] = val


def const(i,j):
    return 1.0

#three point digital filter
def digi3(i,j):
    digi3a = np.array([[ 4., 2.],
                       [ 2., 1.]])
    i2 = abs(i)
    j2 = abs(j)
    if (i2 > 1) or (j2 > 1):
        return 0.0
    else:
        return (0.25/3.0)*digi3a[i2,j2]

#five point digital filter
def digi5(i,j):
    digi5a = np.array([[ 6., 4., 1],
                       [ 4., 0., 0],
                       [ 1., 0., 0]])
    i2 = abs(i)
    j2 = abs(j)
    if (i2 > 2) or (j2 > 2):
        return 0.0
    else:
        return (1.0/16.0)*digi5a[i2,j2]


def init_center(kernel, dfilt=const):
    nx, ny = np.shape(kernel)
    knx = nx/3
    kny = ny/3

    khnx1 = int( floor(1.0*knx/2.0) )
    khny1 = int( floor(1.0*kny/2.0) )
    khnx2 = int(  ceil(1.0*knx/2.0) )
    khny2 = int(  ceil(1.0*kny/2.0) )
    print("size = {} {} ({} {})".format(khnx1, khny1, nx, ny))

    #positive side
    for i in range(0,khnx2+1):
        for j in range(0,khny1+1):
            kernel[i,j] = dfilt(i,j)

    #corner
    for i in range(-1, -khnx1, -1):
        for j in range(-1, -khny1, -1):
            i2 = nx+i
            j2 = ny+j
            #print("  ij=({},{}) i2j2=({},{})".format(i,j,i2,j2))
            kernel[i2,j2] = dfilt(i,j)

    #top
    for i in range(0,khnx2+1):
        for j in range(-1, -khny1, -1):
            i2 = i
            j2 = ny+j
            kernel[i2,j2] = dfilt(i,j)

    #bottom
    for i in range(-1, -khnx1, -1):
        for j in range(0,khny2+1):
            i2 = nx+i
            j2 = j
            kernel[i2,j2] = dfilt(i,j)



class Filters(unittest.TestCase):



    def test_init(self):
        """ write to kernel & image variables and read back"""

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
                self.assertAlmostEqual(kernel[i,j], k2[i,j] , places=4)
                self.assertAlmostEqual( image[i,j],img2[i,j], places=4)
                

    def test_kernel_init(self):
        """ TODO: transform with unitary transformation and get original back"""

        NxMesh = 3
        NyMesh = 3

        #internal filter size is 3x Nx/y/zMesh
        NxF = NxMesh*3
        NyF = NyMesh*3

        flt = pypic.Filter(NxMesh, NyMesh)

        flt.init_kernel()
        kernel = np.zeros((NxF, NyF))
        kernel[0,0] = 1.0
        flt.set_kernel( flatten(kernel) )

        #image  = np.random.rand(NxF, NyF)
        #flt.set_image(  flatten(image) )

        k2   = reshape( flt.get_kernel(), NxF, NyF)
        #img2 = reshape( flt.get_image() , NxF, NyF)
                
        #print()
        #print(k2)

        flt.fft_kernel()
        k3   = reshape( flt.get_kernel(), NxF, NyF)

        #print()
        #print(k3)



    def test_fft_backandforth(self):
        """FFT transform image forward and then backward to see if we get the same result back
         NOTE: floating-point conversion if not exact so we need some error tolerance"""

        NxMesh = 3
        NyMesh = 4 #make Nx != Ny to catch additional index bugs

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
                self.assertAlmostEqual(img1[i,j], img3[i,j], places=4) 



    def skip_smearing_test0(self):
        """ put Gaussian filter in, convolve, and compare to scipy.conv2d"""

        plt.fig = plt.figure(1, figsize=(4,6))
        plt.rc('font', family='serif', size=12)
        plt.rc('xtick')
        plt.rc('ytick')
        
        gs = plt.GridSpec(4, 2)
        
        axs = []
        for ai in range(8):
            axs.append( plt.subplot(gs[ai]) )

        NxMesh = 6
        NyMesh = 6

        #internal filter size is 3x Nx/y/zMesh
        NxF = NxMesh*3
        NyF = NyMesh*3

        vmin = 0.0
        vmax = 1.0

        flt = pypic.Filter(NxMesh, NyMesh)

        #flt.init_kernel()
        #flt.init_gaussian_kernel(5.0)

        #f0 = 2.0/NxMesh/2.0/np.pi
        #f0 = 0.1
        #flt.init_sinc_kernel(f0, f0)

        #flt.init_lowpass_fft_kernel(0)

        kernel = np.zeros((NxF, NyF))
        #kernel[ 0, 0] = 1.0
        #kernel[ 1, 0] = 1.0
        #kernel[ 0, 1] = 1.0
        #kernel[ 1, 1] = 1.0
        #kernel[ 2, 0] = 1.0
        #kernel[ 2, 1] = 1.0
        #kernel[ 1, 2] = 1.0
        #kernel[ 2, 0] = 1.0
        #kernel[ 2, 1] = 1.0
        #kernel[ 2, 2] = 1.0


        #kernel[-1,-1] = 1.0
        #kernel[ 0,-1] = 1.0
        #kernel[ 1,-1] = 1.0
        #kernel[-1, 0] = 1.0
        #kernel[-1, 1] = 1.0

        #kernel[ 2,-1] = 1.0
        #kernel[-1, 2] = 1.0

        #for i in range(0, NxMesh):
        #    for j in range(0,NyMesh):
        #        kernel[i,j] = 1

        print("initing kernel.....")
        #init_center(kernel)
        init_center(kernel, digi3)
        #init_center(kernel, digi5)
        print(kernel)

        #swap kernel
        #(-1)^(i + j + ...)
        #for i in range(NxF):
        #    for j in range(NyF):
        #        kernel[i,j] *= (-1.0)**(i+j)
        #print()
        #print(kernel)
        flt.set_kernel( flatten(kernel) )

        #kernel = np.zeros((NxF, NyF))
        #kernel = sinc2d(kernel, X, Y)
        #flt.set_kernel( flatten(kernel) )

        ##################################################
        image = np.zeros((NxF, NyF))

        #centered on middle
        #image[NxMesh:NxMesh+NxMesh, NyMesh:NyMesh+NyMesh]  = np.random.rand(NxMesh, NyMesh)

        #centered on corner
        #image[0:NxMesh, 0:NyMesh]  = np.random.rand(NxMesh, NyMesh)

        # fill image
        if False:
            for i in range(-2, 1, 1):
                for j in range(-2, 1, 1):
                    iff = fftshift1d(i, NxF)
                    jff = fftshift1d(j, NyF)
                    #print("i={} j={} i2={} j2={}".format(i,j,iff,jff))
                    image[iff,jff] = np.random.rand(1)

        if True:
            for i in range(0, NxMesh):
                for j in range(0, NyMesh):
                    #image[i,j] = np.random.rand(1)
                    image[i+NxMesh,j+NyMesh] = np.random.rand(1)

        flt.set_image(  flatten(image) )

        ##################################################

        ker = reshape( flt.get_kernel(), NxF, NyF)
        img = reshape( flt.get_image( ), NxF, NyF)

        print()
        print(ker)
        print(ker.min(), ker.max())
        print(img.min(), img.max())
                
        axs[0].imshow(ker, vmin=vmin, vmax=vmax)
        axs[1].imshow(img, vmin=vmin, vmax=vmax)

        flt.fft_kernel()
        flt.fft_image_forward()
        ker2 = reshape( flt.get_kernel(), NxF, NyF)
        img2 = reshape( flt.get_image() , NxF, NyF)

        print()
        print(ker2)

        print(ker2.min(), ker2.max())
        print(img2.min(), img2.max())

        axs[2].imshow(fftshift(ker2) )#, vmin=vmin, vmax=vmax)
        axs[3].imshow(fftshift(img2) )#, vmin=vmin, vmax=vmax)
        print(fftshift(ker2))
        print(fftshift(img2))

        #apply kernel
        flt.apply_kernel()
        flt.apply_kernel()
        flt.fft_image_backward()
        ker3 = reshape( flt.get_kernel(), NxF, NyF)
        img3 = reshape( flt.get_image() , NxF, NyF)

        print(ker3.min(), ker3.max())
        print(img3.min(), img3.max())

        axs[4].imshow(ker3, vmin=vmin, vmax=vmax)
        axs[5].imshow(img3, vmin=vmin, vmax=vmax)
        #axs[5].imshow(fftshift(img3), vmin=vmin, vmax=vmax)
        axs[5].imshow(flip(img3), vmin=vmin, vmax=vmax)


        #python conv2d for comparison
        #pimg = img[
        pyimg = convolve2d(img, ker, mode='same', boundary='wrap')
        print(pyimg.min(), pyimg.max())
        axs[6].imshow(fftshift(pyimg), vmin=vmin, vmax=vmax)


        pyimg2 = convolve(img, ker, mode='same' )
        pyimg2 /= pyimg2.max()
        print(pyimg2.min(), pyimg2.max())
        #axs[7].imshow(pyimg2, vmin=vmin, vmax=vmax)

        err = fftshift(pyimg) - flip(img3)
        print(err)
        print("min = {}".format(err.min()))
        print("max = {}".format(err.max()))

        axs[7].imshow(err, vmin=-1.0, vmax=1.0)

        #plt.savefig("filter.png")


    def test_smearing2(self):

        #plt.fig = plt.figure(1, figsize=(4,6))
        #plt.rc('font', family='serif', size=12)
        #plt.rc('xtick')
        #plt.rc('ytick')
        #gs = plt.GridSpec(4, 2)
        #
        #axs = []
        #for ai in range(8):
        #    axs.append( plt.subplot(gs[ai]) )

        NxMesh = 50
        NyMesh = 50

        #internal filter size is 3x Nx/y/zMesh
        NxF = NxMesh*3
        NyF = NyMesh*3

        vmin = 0.0
        vmax = 1.0

        flt = pypic.Filter(NxMesh, NyMesh)

        ###################################################
        # init kernel
        #print("initing kernel.....")
        flt.init_kernel()
        flt.init_3point(1) 

        #kernel = reshape( flt.get_kernel(), NxF, NyF)
        #print(kernel)

        ##################################################
        # create image

        image = np.zeros((NxF, NyF))
        for i in range(0, NxMesh):
            for j in range(0, NyMesh):
                image[i+NxMesh,j+NyMesh] = np.random.rand(1)

        flt.set_image(  flatten(image) )

        ##################################################
        # visualize

        ker = reshape( flt.get_kernel(), NxF, NyF)
        img = reshape( flt.get_image( ), NxF, NyF)
        #print(ker)

        #axs[0].imshow(ker )#, vmin=vmin, vmax=vmax)
        #axs[1].imshow(img, vmin=vmin, vmax=vmax)

        ##################################################
        # fft

        flt.fft_kernel()
        flt.fft_image_forward()
        ker2 = reshape( flt.get_kernel(), NxF, NyF)
        img2 = reshape( flt.get_image() , NxF, NyF)

        #print()
        #print(ker2)

        #print(ker2.min(), ker2.max())
        #print(img2.min(), img2.max())

        #axs[2].imshow(fftshift(ker2) )#, vmin=vmin, vmax=vmax)
        #axs[3].imshow(fftshift(img2) )#, vmin=vmin, vmax=vmax)

        ##################################################
        # apply kernel

        flt.apply_kernel()
        flt.fft_image_backward()
        ker3 = reshape( flt.get_kernel(), NxF, NyF)
        img3 = reshape( flt.get_image() , NxF, NyF)

        #print(ker3.min(), ker3.max())
        #print(img3.min(), img3.max())

        #axs[4].imshow(ker3, vmin=vmin, vmax=vmax)
        #axs[5].imshow(img3, vmin=vmin, vmax=vmax)

        cimg_fft = img3[NxMesh:2*NxMesh, NyMesh:2*NyMesh]
        #axs[5].imshow(cimg_fft, vmin=vmin, vmax=vmax)

        ################################################### 
        # digital filtering for comparison
        flt.set_image( flatten(image) ) # set original image back
        img4 = reshape( flt.get_image( ), NxF, NyF)
        #axs[6].imshow(img4, vmin=vmin, vmax=vmax)

        for i in range(1):
            flt.direct_convolve_3point() # direct convolve

        img5 = reshape( flt.get_image( ), NxF, NyF)

        cimg_dc = img5[NxMesh:2*NxMesh, NyMesh:2*NyMesh]
        #axs[6].imshow(cimg_dc, vmin=vmin, vmax=vmax)
        

        ################################################### 
        # compute error between the two methods
        err = cimg_dc - cimg_fft
        #print( err )
        #axs[7].imshow( err, vmin=-1.0, vmax=1.0)

        #print("max img orig: {}".format(img4.max() ))
        #print("max img fft : {}".format(img3.max() ))
        #print("max img conv: {}".format(img5.max() ))
        #print("max error   : {}".format( np.abs(err).max() ))

        #plt.savefig("filter2.png")

        for i in range(NxMesh):
            for j in range(NyMesh):
                self.assertAlmostEqual(cimg_dc[i,j], cimg_fft[i,j], places=4) 

    def skip_test_gaussian(self):

        plt.fig = plt.figure(1, figsize=(4,6))
        plt.rc('font', family='serif', size=12)
        plt.rc('xtick')
        plt.rc('ytick')
        gs = plt.GridSpec(4, 2)
        
        axs = []
        for ai in range(8):
            axs.append( plt.subplot(gs[ai]) )

        NxMesh = 100
        NyMesh = 100

        #internal filter size is 3x Nx/y/zMesh
        NxF = NxMesh*3
        NyF = NyMesh*3

        vmin = 0.0
        vmax = 1.0

        flt = pypic.Filter(NxMesh, NyMesh)

        ###################################################
        # init kernel and create 3-point image
        flt.init_kernel()
        flt.init_3point(1) 

        image = np.zeros((NxF, NyF))
        for i in range(0, NxMesh):
            for j in range(0, NyMesh):
                image[i+NxMesh,j+NyMesh] = np.random.rand(1)
        flt.set_image(  flatten(image) )

        img = reshape( flt.get_image( ), NxF, NyF)
        axs[0].imshow(img[NxMesh:2*NxMesh, NyMesh:2*NyMesh], vmin=vmin, vmax=vmax)

        flt.fft_kernel()
        flt.fft_image_forward()
        flt.apply_kernel()
        flt.fft_image_backward()

        img2 = reshape( flt.get_image() , NxF, NyF)
        #axs[1].imshow(img2, vmin=vmin, vmax=vmax)

        cimg1 = img2[NxMesh:2*NxMesh, NyMesh:2*NyMesh]
        axs[1].imshow(cimg1, vmin=vmin, vmax=vmax)




        ##################################################
        # Gaussian (in FFT space) for comparison

        flt2 = pypic.Filter(NxMesh, NyMesh)
        flt2.set_image(  flatten(image) )

        sigx = 6.5
        sigy = 6.5
        #flt2.init_gaussian_kernel(sigx, sigy)
        flt2.init_lowpass_fft_kernel(140)

        #plot kernel in Fourier space
        ker2 = reshape( flt.get_kernel(), NxF, NyF)
        img3 = reshape( flt2.get_image() , NxF, NyF)

        cker1 = reshape( flt.get_kernel(), NxF, NyF)

        axs[2].imshow(img3[NxMesh:2*NxMesh, NyMesh:2*NyMesh])#, vmin=vmin, vmax=vmax)
        axs[3].imshow(fftshift(ker2) )#, vmin=vmin, vmax=vmax)


        flt2.fft_image_forward()
        flt2.apply_kernel()
        flt2.fft_image_backward()

        # visualize
        ker4 = reshape( flt2.get_kernel(), NxF, NyF)
        img4 = reshape( flt2.get_image() , NxF, NyF)

        cimg2 = img4[NxMesh:2*NxMesh, NyMesh:2*NyMesh]
        cker2 = fftshift(ker4)
        axs[4].imshow(cimg2, vmin=vmin, vmax=vmax)
        axs[5].imshow(cker2 )#, vmin=vmin, vmax=vmax)



        ################################################### 
        # compute error between the two methods
        ierr = (cimg1 - cimg2)/cimg1
        #print( err )
        axs[6].imshow( ierr, vmin=-1.0, vmax=1.0)

        kerr = (cker1 - cker2)/cker1
        axs[7].imshow( kerr, vmin=-1.0, vmax=1.0)

        #print("max img orig: {}".format(img4.max() ))
        #print("max img fft : {}".format(img3.max() ))
        #print("max img conv: {}".format(img5.max() ))
        #print("max error   : {}".format( np.abs(err).max() ))

        #for i in range(NxMesh):
        #    for j in range(NyMesh):
        #        self.assertAlmostEqual(cimg_dc[i,j], cimg_fft[i,j], places=10) 


        print("max img orig: {}".format(img.max() ))
        print("max img fft : {}".format(img2.max() ))
        print("max img gaus: {}".format(img4.max() ))
        print("max error   : {}".format( np.abs(ierr).max() ))



        #plt.savefig("filter3.png")

    def test_filters_in_action(self):

        #plt.fig = plt.figure(1, figsize=(4,6))
        #plt.rc('font', family='serif', size=12)
        #plt.rc('xtick')
        #plt.rc('ytick')
        #gs = plt.GridSpec(3, 1)
        #
        #axs = []
        #for ai in range(3):
        #    axs.append( plt.subplot(gs[ai]) )


        conf = Conf()
        conf.Nx = 3
        conf.Ny = 3
        conf.Nz = 1
        conf.NxMesh = 10
        conf.NyMesh = 10
        conf.NzMesh = 1

        node = pycorgi.Node2D(conf.Nx, conf.Ny)
        node.setGridLims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)
        loadTiles(node, conf)
        insert_em(node, conf, linear_ramp)
        #inject(node, filler_no_velocity, conf) #injecting plasma particles

        flt = pypic.Filter(conf.NxMesh, conf.NyMesh)
        flt.init_gaussian_kernel(4.0, 4.0)

        flt.get_padded_current( node.getTilePtr(1,1), node)

        img = reshape( flt.get_image( ), conf.NxMesh*3, conf.NyMesh*3)
        #axs[0].imshow(img[conf.NxMesh:2*conf.NxMesh, conf.NyMesh:2*conf.NyMesh]) #, vmin=vmin, vmax=vmax)
        #axs[0].imshow(img, vmin=0.0, vmax=100.0) #, vmin=vmin, vmax=vmax)


        # reference array
        data = np.zeros((conf.Nx*conf.NxMesh, conf.Ny*conf.NyMesh, conf.Nz*conf.NzMesh, 3))
        for cid in node.getTileIds():
            c = node.getTilePtr( cid )
            (i, j) = c.index()

            yee = c.getYee(0)
            for k in range(conf.NyMesh):
                for q in range(conf.NxMesh):
                    for r in range(conf.NzMesh):
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 0] = yee.jx[q,k,r]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 1] = yee.jy[q,k,r]
                        data[ i*conf.NxMesh + q, j*conf.NyMesh + k, 0*conf.NzMesh + r, 2] = yee.jz[q,k,r]

        #img2 = data[:,:,0,0]
        #axs[1].imshow(img2, vmin=0.0, vmax=100.0) #, vmin=vmin, vmax=vmax)

        for i in range(0,3*conf.NxMesh):
            for j in range(0,3*conf.NyMesh):
                self.assertEqual(data[i,j,0,0], img[i,j]) #jx

        flt.fft_image_forward()
        flt.apply_kernel()
        flt.fft_image_backward()

        #img = reshape( flt.get_image( ), conf.NxMesh*3, conf.NyMesh*3)
        #axs[2].imshow(img, vmin=0.0, vmax=100.0) #, vmin=vmin, vmax=vmax)

        #plt.savefig("filter_getter.png")
