import numpy as np
import os, sys
import copy
import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm

#from visualize import *
import mcmc

#set seed so we get predictable results
np.random.seed( 12 )
from mpl_toolkits.mplot3d import Axes3D


#Profiling / Timing
sys.path.insert(0, '../tools')
from timer import Timer



class Params:
    mins = None
    maxs = None
    lens = None


def randab(a, b):
    return a + (b-a)*np.random.rand()


#random 2d location
def randLoc(params):
    x = randab(params.mins[0], params.maxs[0])
    y = randab(params.mins[1], params.maxs[1])
    return (x,y)


#random 3D direction in terms of angles phi and theta
def randUnitSphere():
    vphi = randab( 0.0, 2.0*np.pi ) #azimuth
    vthe = randab(-np.pi/2.0, np.pi/2.0 ) #latitude
    return (vphi, vthe)

def sph2cart(vr, vphi, vthe):
    vx = vr * np.sin(vthe) * np.cos(vphi)
    vy = vr * np.sin(vthe) * np.sin(vphi)
    vz = vr * np.cos(vthe)
    return vx, vy, vz

def pol2cart(vr, vphi, vthe):
    vx = vr * np.cos(vphi)
    vy = vr * np.sin(vphi)
    return vx, vy, 0.0

def randVel(vabs):
    (vphi, vthe) = randUnitSphere()
    vx, vy, vz = sph2cart(vabs, vphi, vthe)
    #vx, vy, vz = pol2cart(vabs, vphi, vthe)
    return vx, vy, vz


def unitVecX(vec):
    return np.array([ 0.0, vec[0] ])
def unitVecY(vec):
    return np.array([ 0.0, vec[1] ])
def unitVecZ(vec):
    return np.array([ 0.0, vec[2] ])

def norm(vec):
    return np.sqrt( vec[0]**2 + vec[1]**2 + vec[2]**2 )

#unit vector into direction a x b
def uCross(vecA, vecB):
    vecC = np.cross(vecA, vecB)
    return vecC / norm(vecC)

def printVec(vec, name):
    print " {}: ({}, {}, {}) |v| = {} ".format(name, vec[0], vec[1], vec[2], norm(vec))


########################################

# Monte Carlo Compton scattering according to Sobol 1977
def comptonScatter(e, p, ax, plot=True):

    # notes about units
    # p.hv = hv/m_e c^2
    # e.v = v/c 

    beta = np.array([ e.vx(), e.vy(), e.vz() ]) / e.beta() #electron unit vector
    omeg = np.array([ p.vx(), p.vy(), p.vz() ])         #photon unit vector
    theta = np.arccos( np.dot(beta, omeg) ) #angle between electron and photon


    #printVec(beta, "beta")
    #printVec(omeg, "omeg")
    #print "thet:", theta


    #rotate to scattering plane (i,j,k)
    kvec = np.array([ 0.0, -1.0, 0.0] ) #this points along r, here we select it arbitrarily
    #kvec = omeg #another choice so we are in direction of photon
    jvec = uCross(beta, omeg)
    ivec = uCross(kvec, jvec)
    M = np.array([ ivec, jvec, kvec ])

    #printVec(ivec, "i")
    #printVec(jvec, "j")
    #printVec(kvec, "k")


    ################################################## 
    #print M
    #print np.linalg.inv(M)


    ################################################## 

    #unit vector of electron in scattering coordinates
    v0  = np.matmul(M, np.array([ e.vx(), e.vy(), e.vz()] )) #rotate electron velocity to scattering plane (i,k)
    mu  = v0[0]*np.sin(theta) +  v0[2]*np.cos(theta)
    rho = np.sqrt( v0[0]**2 + v0[1]**2 )

    #printVec(v0, "v0")
    
    #Compton parameters
    x = 2.0 * p.hv() * e.gamma() * (1.0 - mu*e.beta() )
    y = x/2.0

    #print "Compton x: {}".format(x)

    #additional scattering angles (v0, w0, t0) define a frame of reference
    w0 = np.array([ v0[1], -v0[0], 0])/rho
    t0 = np.array([ v0[0]*v0[2], v0[1]*v0[2], -rho**2])/rho

    #printVec(w0, "w0")
    #printVec(t0, "t0")

    #ax.plot( unitVecX(v0), unitVecY(v0), unitVecZ(v0), label='v0' )
    #ax.plot( unitVecX(w0), unitVecY(w0), unitVecZ(w0), label='w0' )
    #ax.plot( unitVecX(t0), unitVecY(t0), unitVecZ(t0), label='t0' )


    #scatter
    done = False
    OmegaOmegap = 0.0
    while not(done):
        z1 = np.random.rand()
        z2 = np.random.rand()
        z3 = np.random.rand()

        mup  = (e.beta() + 2.0*z1 - 1.0)/(1.0 + e.beta()*(2.0*z1 - 1.0))
        phip = 2.0*np.pi*z2

        OmegaOmegap = mu*mup - np.sqrt( 1.0 - mup**2) * ( rho*np.sin(phip)*np.cos(theta) 
                - (1.0/rho)*(v0[1]*np.cos(phip) + v0[0]*v0[2]*np.sin(phip))*np.sin(theta))

        yp = y / (1.0 + p.hv()*(1.0 - OmegaOmegap))/(e.gamma() * (1.0 - mup*e.beta())) 

        Y = yp/y + (yp/y)**3 + (yp/y)**2 *( (1.0/yp - 1.0/y)**2 - 2.0*( 1.0/yp - 1.0/y) )

        if Y > 2.0*z3:
            done = True
    #now we have scattered successfully

    #new energy
    hvp = yp / ( e.gamma()*(1.0 - mup*e.beta()) ) 
    #print " energy shift: {}".format( hvp/p.hv() )

    #new direction
    Omegap_ijk = mup*v0 + np.sqrt(1.0-mup**2)*( w0*np.cos(phip) + t0*np.sin(phip) )
    Omegap = np.dot( np.linalg.inv(M), Omegap_ijk ) 
    Omegap = Omegap / np.sqrt( Omegap[0]**2 + Omegap[1]**2 + Omegap[2]**2 ) #normalize

    if plot:
        ax.plot( unitVecX(Omegap), unitVecY(Omegap), unitVecZ(Omegap), alpha=0.2, color='red')#, label='Omegap' )

    phs = mcmc.photon( hvp, Omegap[0], Omegap[1], Omegap[2] )


    return e, phs



def visualize(ax, slab, params, lap):
    ax.cla()

    # Normalize 3d plot to make aspect ratio unity
    #MAX = 1.0
    #for direction in (-1, 1):
    #    for point in np.diag(direction * MAX * np.array([1,1,1])):
    #        axs[0].plot([point[0]], [point[1]], [point[2]], 'w')

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim( 0.0, 1.0)

    xs = slab.xloc
    ys = slab.yloc
    zs = slab.zloc

    ax.scatter(xs, ys, zs, 
            color='black',
            marker='.'
            )
            

    slap = str(lap).rjust(4, '0')
    fname = 'slab_'+slap+'.png'
    plt.savefig(fname)

    return 

def angleHist(ax, slab):
    ax.cla()

    xmin = 0.0
    xmax = 1.0
    #ax.set_xscale('log')
    #ax.set_yscale('log')

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0.0, 1.1)
    ax.minorticks_on()
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel(r'I')
    
    vz = slab.overflow.vecZ() #z component of velocity (|v| = 1.0, i.e., normalized)

    #theta = np.cos(vz)
    #mu = np.acos(theta)
    mu = vz

    hist, edges = np.histogram(mu, np.linspace(xmin, xmax, 20))
    hist = 1.0 * hist / hist.max() #normalize
    ax.plot(edges[:-1], hist )



def bbodyrad(x, kT):
    return (15.0/np.pi**4) * ((x/kT)**3)/(np.exp(x/kT) - 1)


def energyHist(ax, slab):

    #prepare axis
    ax.cla()

    xmin = 0.5
    xmax = 100.0
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0.2, 1.1)
    ax.minorticks_on()

    #histogram
    xs = slab.overflow.vecE()
    hist, edges = np.histogram(xs, np.logspace(np.log10(xmin), np.log10(xmax), 50))
    hist = 1.0 * hist / hist.max() #normalize
    ax.plot(edges[:-1], hist )

    #test against real Planck
    kT = 1.0
    xx = np.logspace(np.log10(xmin), np.log10(xmax), 50)
    yy = bbodyrad(xx, kT)
    yy = yy / yy.max()
    ax.plot(xx, yy, "r-")

    ax.plot([2.95*kT, 2.95*kT], [0.0, 1.0], "g--")




def electronHist(ax, slab):

    #prepare axis
    ax.cla()


    if True:
        xmin = 0.8
        xmax = 10.0
        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0.00001, 1.1)
        ax.minorticks_on()
        ax.set_xlabel(r'$\gamma$')
        ax.set_ylabel(r'$f(\gamma)$')


        #draw samples from electron distribution
        zs = np.zeros(100000)
        for i in range(len(zs)):
            el = slab.sampleBoostedMaxw(0.201, [0.0, 0.0, 0.0])
            zs[i] = el.gamma()
            #print el.gamma()


        #histogram
        hist, edges = np.histogram(zs, np.logspace(np.log10(xmin), np.log10(xmax), 50))
        #hist, edges = np.histogram(zs, np.linspace(xmin, xmax, 50))
        hist = 1.0 * hist / hist.max() #normalize
        ax.plot(edges[:-1], hist )


        if False:
            zs = np.zeros(100000)
            for i in range(len(zs)):
                el = slab.sampleBoostedMaxw(0.20, [0.0, 0.0, 0.0])
                zs[i] = el.gamma()

            #histogram
            hist, edges = np.histogram(zs, np.logspace(np.log10(xmin), np.log10(xmax), 50))
            hist = 1.0 * hist / hist.max() #normalize
            ax.plot(edges[:-1], hist, "r-")


    # beta
    if False:
        xmin = 0.0
        xmax = 1.0
        #ax.set_xscale('log')
        #ax.set_yscale('log')

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0.0, 1.1)
        ax.minorticks_on()
        ax.set_xlabel(r'$\beta$')
        ax.set_ylabel(r'$f(\beta)$')


        #draw samples from electron distribution
        zs = np.zeros(100000)
        for i in range(len(zs)):
            el = slab.sampleBoostedMaxw(0.101, [0.0, 0.0, 0.0])
            #zs[i] = el.beta()
            zs[i] = el.vx()
            #zs[i] = np.sqrt(el.vx()**2 + el.vy()**2 + el.vz()**2)

        #histogram
        #hist, edges = np.histogram(zs, np.logspace(np.log10(xmin), np.log10(xmax), 50))
        hist, edges = np.histogram(zs, np.linspace(xmin, xmax, 50))
        hist = 1.0 * hist / hist.max() #normalize
        ax.plot(edges[:-1], hist )

        print "rel std: {}".format(np.std(zs))

        zs = np.zeros(100000)
        for i in range(len(zs)):
            el = slab.sampleBoostedMaxw(0.10, [0.0, 0.0, 0.0])
            #zs[i] = el.beta()
            zs[i] = el.vx()
            #zs[i] = np.sqrt(el.vx()**2 + el.vy()**2 + el.vz()**2)


        #hist, edges = np.histogram(zs, np.logspace(np.log10(xmin), np.log10(xmax), 50))
        hist, edges = np.histogram(zs, np.linspace(xmin, xmax, 50))
        hist = 1.0 * hist / hist.max() #normalize
        ax.plot(edges[:-1], hist, "r-")

        print "non-rel std: {}".format(np.std(zs))

        kT = 0.1

        xx = np.linspace(0.0, 1.0, 100)
        yy = np.exp(-xx**2/(kT**2))

        ax.plot(xx, yy, "k-")




if __name__ == "__main__":


    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')

    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(4,7))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(3, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0], projection='3d') )
    axs.append( plt.subplot(gs[1]) )
    axs.append( plt.subplot(gs[2]) )


    ################################################## 
    # set-up grid
    # xy
    params = Params()
    params.mins = [-1.0,-1.0, 0.0 ]
    params.maxs = [ 1.0, 1.0, 1.0 ]
    params.lens = [ 1.0, 1.0, 1.0 ]


    timer = Timer(["bucket", "step"])

    ################################################### 
    # Create bucket
    bucket = mcmc.photonBucket()
    print "created bucket ({})".format( bucket.size() )

    #pour isotropic photons to the bucket for testing
    timer.start("bucket")
    for i in range(10):
        (vx, vy, vz) = randVel(1.0) #direction on unit sphere
        E = 0.01 # E = h\nu 

        ph = mcmc.photon( E, vx, vy, vz )
        bucket.push_back( ph )
    print "loaded bucket with {} photons".format( bucket.size() )
    timer.stats("bucket")


    ##################################################
    #Lets create a simulation slab
    slab = mcmc.Slab( bucket )
    slab.set_dimensions(
            params.mins[0], params.maxs[0],
            params.mins[1], params.maxs[1],
            params.mins[2], params.maxs[2]
                       )
    slab.set_numberDensity(0.5)

    slab.floor() # put everything to the bottom of the slab

    #flux from the bottom
    flux = 100.0 #TODO units

    timer.start("step")
    #slab.inject(1000.0)

    for lap in range(50):
        print "----lap: {}".format(lap)

        slab.inject(flux)
        slab.push()
        slab.wrap()
        slab.scrape()
        #slab.scatter(0.21)

        timer.lap("step")

        #angleHist(axs[1], slab)
        #energyHist(axs[2], slab)
        #visualize(axs[0], slab, params, lap)
    print "overflow size: {}".format(slab.overflow.size())

    timer.stop("step")
    timer.stats("step")

    angleHist(axs[1], slab)
    #energyHist(axs[2], slab)
    electronHist(axs[2], slab)
    visualize(axs[0], slab, params, 100)
    


    plt.savefig("slab.png")
