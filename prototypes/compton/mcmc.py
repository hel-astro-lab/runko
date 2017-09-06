import numpy as np
import os, sys
import copy
import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm

#from visualize import *
import mcmc

#set seed so we get predictable results
np.random.seed( 1234 )
from mpl_toolkits.mplot3d import Axes3D


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
    return (vthe, vphi)

def sph2cart(vr, vphi, vthe):
    vx = vr * np.sin(vthe) * np.cos(vphi)
    vy = vr * np.sin(vthe) * np.sin(vphi)
    vz = vr * np.cos(vthe)
    return vx, vy, vz


def randVel(vabs):
    (vphi, vthe) = randUnitSphere()
    vx, vy, vz = sph2cart(vabs, vphi, vthe)
    return vx, vy, vz


def unitVecX(vec):
    return np.array([ 0.0, vec[0] ])
def unitVecY(vec):
    return np.array([ 0.0, vec[1] ])
def unitVecZ(vec):
    return np.array([ 0.0, vec[2] ])


########################################

# Monte Carlo Compton scattering according to Sobol 1977
def comptonScatter(e, p, ax, plot=True):

    # notes about units
    # p.hv = hv/m_e c^2
    # e.v = v/c = beta

    beta = np.array([ e.vx(), e.vy(), e.vz() ]) / e.v() #electron unit vector
    omeg = np.array([ p.vx(), p.vy(), p.vz() ])         #photon unit vector
    theta = np.arccos( np.dot(beta, omeg) ) #angle between electron and photon

    #rotate to scattering plane (i,j,k)
    kvec = beta
    jvec = np.cross(beta, omeg)
    ivec = np.cross(kvec, jvec)
    M = np.array([ ivec, jvec, kvec ])

    ################################################## 
    #print M
    #print np.linalg.inv(M)

    #MAX = 3
    #for direction in (-1, 1):
    #    for point in np.diag(direction * MAX * np.array([1,1,1])):
    #        ax.plot([point[0]], [point[1]], [point[2]], 'w')
    #ax.plot( unitVecX(beta), unitVecY(beta), unitVecZ(beta), label='beta' )
    #ax.plot( unitVecX(omeg), unitVecY(omeg), unitVecZ(beta), label='omeg' )
    #ax.plot( unitVecX(ivec), unitVecY(ivec), unitVecZ(ivec), label='i' )
    #ax.plot( unitVecX(jvec), unitVecY(jvec), unitVecZ(jvec), label='j' )
    #ax.plot( unitVecX(kvec), unitVecY(kvec), unitVecZ(kvec), label='k' )

    ################################################## 

    #unit vector of electron in scattering plane
    v0  = np.matmul(M, beta)
    mu  = v0[0]*np.sin(theta) +  v0[2]*np.cos(theta)
    rho = np.sqrt( v0[0]**2 + v0[1]**2 )
    
    x = 2.0 * p.hv() * e.gamma() * (1.0 - mu*e.v() )
    y = x/2.0

    #additional scattering angles
    w0 = np.array([ v0[1], -v0[0], 0])/rho
    t0 = np.array([ v0[1]*v0[2], v0[1]*v0[2], -rho**2])/rho

    # draw scatter
    done = False
    OmegaOmegap = 0.0
    while not(done):
        z1 = np.random.rand()
        z2 = np.random.rand()
        z3 = np.random.rand()

        mup  = (e.v() + 2.0*z1 - 1.0)/(1.0 + e.v()*(2.0*z1 - 1.0))
        phip = 2.0*np.pi*z2

        OmegaOmegap = mu*mup - np.sqrt( 1.0 - mup**2) * ( rho*np.sin(phip)*np.cos(theta) 
                - (1.0/rho)*(v0[1]*np.cos(phip) + v0[0]*v0[2]*np.sin(phip))*np.sin(theta))

        yp = y / (1.0 + p.hv()*(1.0 - OmegaOmegap))/(e.gamma() * (1.0 - mup*e.v())) 

        Y = yp/y + (yp/y)**3 + (yp/y)**2 *( (1.0/yp - 1.0/y)**2 - 2.0*( 1.0/yp - 1.0/y) )

        if Y > 2.0*z3:
            done = True
    #now we have scattered successfully

    #new energy
    hvp = yp / ( e.gamma()*(1.0 - mup*e.v()) ) 

    #print " energy shift: {}".format( hvp/p.hv() )


    #new direction
    Omegap_ijk = mup*v0 + np.sqrt(1.0-mup**2)*( w0*np.cos(phip) + t0*np.sin(phip) )
    Omegap = np.dot( np.linalg.inv(M), Omegap_ijk ) 
    Omegap = Omegap / np.sqrt( Omegap[0]**2 + Omegap[1]**2 + Omegap[2]**2 ) #normalize

    if plot:
        ax.plot( unitVecX(Omegap), unitVecY(Omegap), unitVecZ(Omegap) )#, label='Omegap' )

    phs = mcmc.photon( hvp, Omegap[0], Omegap[1], Omegap[2] )

    return e, phs



if __name__ == "__main__":


    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')

    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,4))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 2)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1], projection='3d') )
    


    ################################################## 
    # set-up grid
    # xy
    params = Params()
    params.mins = [ 0.0, 0.0, 0.0 ]
    params.maxs = [ 1.0, 1.0, 1.0 ]
    params.lens = [ 1.0, 1.0, 1.0 ]


    ################################################## 
    bucket = mcmc.photonBucket()
    print "created bucket ({})".format( bucket.size() )

    #pour isotropic photons to the bucket
    for i in range(1000):
        (vx, vy, vz) = randVel(0.8) #direction on unit sphere
        E = 1.0 # E = h\nu 

        ph = mcmc.photon( E, vx, vy, vz )
        bucket.push_back( ph )
    print "loaded bucket with {} photons".format( bucket.size() )


    #isotropic electron with beta = 0.5
    (vx, vy, vz) = randVel(0.5) 
    e = mcmc.electron(1.0, vx, vy, vz)
    print "target electron with beta: {} gamma: {}".format(e.v(), e.gamma() )

    

    print "Compton scatter single interaction..."
    (vx, vy, vz) = randVel(1.0) #direction on unit sphere
    ph = mcmc.photon( 1.0, vx, vy, vz )

    for sc in range(100):
        es, ps = comptonScatter(e, ph, axs[1], plot=True)

    #plt.show()
    #axs[1].legend()

    ################################################## 
    # now Monte Carlo scatter whole bucket

    #axs[0].set_xscale('log')
    #axs[0].set_yscale('log')

    axs[0].set_xlim(0.01, 12.0)
    axs[0].set_ylim(0.1, 1.1)

    #step in time
    xs = np.zeros( bucket.size() )
    for lap in range(5):
        for i in range(bucket.size() ):
            
            #draw random electron
            (vx, vy, vz) = randVel( 0.9  ) 
            e = mcmc.electron(1.0, vx, vy, vz)

            p = bucket.get(i)
            xs[i] = p.hv()

            es, ps = comptonScatter(e, p, axs[1], plot=False)

            bucket.replace(i, ps)

            
    hist, edges = np.histogram(xs, np.linspace(1.1, 2.0, 50))
    hist = hist / np.max(hist)
    axs[0].plot(edges[:-1], hist, "k-")


    #visualize_data(axs[0], mesh, params)
    plt.savefig("compton.png")
