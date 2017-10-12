import numpy as np
import os, sys
import copy
import matplotlib.pyplot as plt
#import palettable as pal
from matplotlib import cm

#from visualize import *
import mcmc

#set seed so we get predictable results
np.random.seed( 12 )
from mpl_toolkits.mplot3d import Axes3D


#Profiling / Timing
sys.path.insert(0, '../tools')
from timer import Timer

sigma_T = 8.0*3.14/3.0    # Thomson cross-section = 8 \pi r_e / 3


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

def scaleVecX(vec):
    return np.array([ 0.0, vec[0] ])
def scaleVecY(vec):
    return np.array([ 0.0, vec[1] ])
def scaleVecZ(vec):
    return np.array([ 0.0, vec[2] ])


def norm(vec):
    return np.sqrt( vec[0]**2 + vec[1]**2 + vec[2]**2 )

#unit vector into direction a x b
def uCross(vecA, vecB):
    vecC = np.cross(vecA, vecB)
    if norm(vecC) == 0.0: return 0.0
    return vecC / norm(vecC)

def printVec(vec, name):
    print " {}: ({}, {}, {}) |v| = {} ".format(name, vec[0], vec[1], vec[2], norm(vec))



########################################

# Monte Carlo Compton scattering according to Sobol 1977
# The reference frame where the scattering is calculated is chosen as:
# k = Omega (direction of incoming photon)
# j = k x beta (incoming electron is in the (i,k) plane)
# i - complementing the right-handed basis
def comptonScatter(e, p, ax, plot):

    # notes about units
    # p.hv = hv/m_e c^2
    # e.v = v/c 

    beta0 = np.array([ e.vx(), e.vy(), e.vz() ]) / e.v() # unit vector in direction of electron motion in lab frame
    omega = np.array([ p.vx(), p.vy(), p.vz() ])         # unit vector in photon direction in lab frame

    #choose scattering frame (i,j,k)
    kvec = omega # k vector along photon direction omega
    jvec = uCross(kvec, beta0)
    ivec = uCross(jvec, kvec)
    M = np.array([ ivec, jvec, kvec ])  # transformation matrix between lab frame and ijk frame

    
#
#   To add exception when beta || omega
#
#    if norm(jvec) == 0.0:
#        jvec = uCross(kvec,beta)
#        ivec = uCross(jvec, kvec)
    cosalpha = np.dot( kvec, beta0 )    # cosine of the angle between electron and k-vector
    sinalpha = np.sqrt( 1.0 - cosalpha**2 ) # sine of the same angle
  
    ################################################## 

    mu = cosalpha  # in ijk frame angle between e and ph equals to the angle between e and k vector
        
    y = p.hv() * e.gamma() * (1.0 - mu*e.v() ) # initial photon and electron 4-product
    
    #scatter
    done = False
    OmegaOmegap = 0.0
    while not(done):
        z1 = np.random.rand()
        z2 = np.random.rand()
        z3 = np.random.rand()
        
        # draw new possible angles
        mup  = (e.v() + 2.0*z1 - 1.0)/(1.0 + e.v()*(2.0*z1 - 1.0))  # cos(alpha') = k \dot Omega'
        phip = 2.0*np.pi*z2               # azimuthal angle calculated from (-j)
        sinalphap = np.sqrt( 1.0 - mup**2 ) #

        OmegaOmegap = mu * mup - sinalphap *  np.sin(phip) * sinalpha  # angle between incoming and outgoing photons
        
        yp = y / (1.0 + p.hv() * (1.0 - OmegaOmegap) / ( e.gamma() * (1.0 - mup*e.v()) ))

        YY = yp/y + (yp/y)**3 + (yp/y)**2 *( (1.0/yp - 1.0/y)**2 - 2.0*( 1.0/yp - 1.0/y) )
 
        if YY > 2.0*z3:
            done = True
    #now we have scattered successfully

    #new energy
    hvp = yp / ( e.gamma()*(1.0 - mup*e.v()) ) 
    
    #hvp2 = p.hv() * (1. - e.v() * mu) / (1. - e.v() * mup + p.hv() * (1. - OmegaOmegap)/ e.gamma()) # energy test
    
    #print"hv = {}, hvp2 ={}".format(hvp,hvp2)   # compare new photon energy calculated different ways
    #print"hvp*p.hv*(1-OOp) = {}, y-yp = {}".format(hvp*p.hv()*(1.-OmegaOmegap),y-yp) # check if the quantities are conserved
    
    
    #new direction in ijk coordinate system
    Omegap_ijk = np.array( [mup*sinalpha + sinalphap * np.sin(phip) * mu, -sinalphap * np.cos(phip),
                            mup * cosalpha - sinalphap * np.sin(phip) * sinalpha] )
    
    Omegap = np.dot( np.linalg.inv(M), Omegap_ijk ) # transferring back to lab system

    phs = mcmc.photon( hvp, Omegap[0], Omegap[1], Omegap[2] )

    # scattered electron parameters
    gammaes = e.gamma() + p.hv() - hvp		
    vxes = ( e.gamma() * e.vx() + p.hv() * omega[0] - hvp * Omegap[0] ) / gammaes
    vyes = ( e.gamma() * e.vy() + p.hv() * omega[1] - hvp * Omegap[1] ) / gammaes
    vzes = ( e.gamma() * e.vz() + p.hv() * omega[2] - hvp * Omegap[2] ) / gammaes
    es = mcmc.electron( gammaes, vxes, vyes, vzes )
    ves= np.array([ vxes, vyes, vzes ])
    
    if plot:
        ax.plot( scaleVecX(hvp*Omegap), scaleVecY(hvp*Omegap), scaleVecZ(hvp*Omegap), alpha=0.2, linestyle='solid', color='red')#, label='Omegap' )
        ax.plot( scaleVecX(ves), scaleVecY(ves), scaleVecZ(ves), alpha=0.2, linestyle='dashed', color='blue')#, label='Omegap' )

    return es, phs	### ! returns scattered electron and photon
#   return e, phs # return scattered photon and initial electron



### !!! Being tested at the moment !!!

## Monte-Carlo Compton scattering for bucket of electrons and photons
## Choose photon, compute whether the scattering occurs at all using maximal (Thomson) cross-section, 
## choose electron to interact with, compute real cross-section, do scattering (by calling comptonScatter)
def lp_ComptonScatter(bucket_el, bucket_ph, axs, deltat, V):
    
    global sigma_T
    
    hv_sc = np.zeros(bucket_ph.size())
    
    
    P_it_max = 2.0 * sigma_T / V
    #print"P_max = {}, bucket size = {}".format(P_it_max,bucket_ph.size())
    
    
    
    for i in range( bucket_ph.size() ):
        
        
        p = bucket_ph.get(i)
        hv_sc[i] = p.hv()

        z1 = np.random.rand()
        t_it = -np.log( z1 ) /  P_it_max / bucket_el.size()

        if deltat < t_it : continue   # if timestep is less than mean free path, leave this photon and go to next one

        j = int(np.floor( np.random.rand() * bucket_el.size() ))
        #print"j = {}".format(j)
        e = bucket_el.get(j)
        
        beta = np.array([ e.vx(), e.vy(), e.vz() ])  # beta in lab frame
        beta0 = np.array([ e.vx(), e.vy(), e.vz() ]) / e.v() # unit vector in direction of electron motion in lab frame
        omega = np.array([ p.vx(), p.vy(), p.vz() ])         # unit vector in photon direction in lab frame
        mu = np.dot(beta0, omega)    # cosine of angle between incident photon and electron

        #calculate real cross-section and real probability of interaction

        xi = e.gamma() * p.hv() * (1.0 - e.v() * mu) # product of initial electron and photon 4-momenta
        
        if xi < 0.01: # Approximate Taylor expansion for Compton total cross-section if xi<0.01, error about 5e-6
            s0 = 1.0 - 2.0 * xi + 5.2 * xi**2 - 9.1 * xi**3 
                # + 1144.0 * xi**4 / 35.0  - 544.0 * xi**5 / 7. + 1892.0 * xi**6 / 21.0   # higher-order terms, not needed if xi<0.01 
        else:    # Exact formula for Klein - Nishina cross section in units of sigma_T
            s0 = 3.0 / 8.0 / xi**2 * ( 4.0 + (xi - 2.0 - 2.0 / xi) * np.log(1.0 + 2.0 * xi) + 
                                       2.0 * xi**2 * (1.0 + xi) / (1.0 + 2.0 * xi)**2 )
        
        vrel = 1.0 - e.v() * mu  # relative direction ("v_rel" in Stern et al. 1995)
        P_it_real = s0 * vrel / 2.0
        
        z2 = np.random.rand()
        if z2 < P_it_real: 
            es, ps = comptonScatter(e, p, axs[1], False)
            bucket_el.replace(j, es)
            bucket_ph.replace(i, ps)
            hv_sc[i] = ps.hv()
            #print"# {} scattered successfully, new energy = {}".format(i, hv_sc[i])
        else: 
            continue
                
        
    hist, edges = np.histogram(hv_sc, np.logspace(np.log10(xmin), np.log10(xmax), 20))
    #print"hist = {}".format(hist)
    hist = 1.0 * hist / hist.max()
    axs[0].plot(edges[:-1], hist )#, "k-")

        
    return bucket_el, bucket_ph
    
    

def visualize(ax, slab, params):

    xs = slab.xloc
    ys = slab.yloc
    zs = slab.zloc

    ax.plot(xs, ys, zs, 
            color='black',
            marker='.'
            )
            

    return 



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
    #axs.append( plt.subplot(gs[0], projection='3d') )
    plt.xlabel('x')
    plt.ylabel('y')
#    plt.zlabel('z')
    

    ################################################## 
    # Normalize 3d plot to make aspect ratio unity
    MAX = 1.0
    for direction in (-1, 1):
        for point in np.diag(direction * MAX * np.array([1,1,1])):
            #axs[1].plot([point[0]], [point[1]], [point[2]], 'w')
            #axs[2].plot([point[0]], [point[1]], [point[2]], 'w')

            axs[0].plot([point[0]], [point[1]], [point[2]], 'w')

    ################################################## 
    # set-up grid
    # xy
    params = Params()
    params.mins = [ 0.0, 0.0, 0.0 ]
    params.maxs = [ 1.0, 1.0, 1.0 ]
    params.lens = [ 1.0, 1.0, 1.0 ]

#    sys.exit()
    ################################################### 
    #test scattering for single case
    print "Compton scatter single interaction..."

    verand=0.8
    (vx, vy, vz) = [0.0, 0.95, 0.0] #randVel(verand) #(0.1, 0.0, 0.0) #randVel(0.5) 
    e = mcmc.electron(1.0, vx, vy, vz)
    print "target electron with beta: {} gamma: {}".format(e.v(), e.gamma() )

    (vx, vy, vz) = [1.0, 0.0, 0.0] #randVel(1.0) #direction on unit sphere (-1.0, 0.0, 0.0) 
    ph = mcmc.photon( 1.0e-4, vx, vy, vz )


    #visualize starting point of collision
    beta = np.array([ e.vx(), e.vy(), e.vz() ]) / e.v() #electron unit vector
    omeg = np.array([ ph.vx(), ph.vy(), ph.vz() ])         #photon unit vector
#    axs[1].plot( unitVecX(-beta), unitVecY(-beta), unitVecZ(-beta), linestyle='dashed', color='black', linewidth=2.0)
    axs[1].plot( scaleVecX(-beta), scaleVecY(-beta), scaleVecZ(-beta), linestyle='dashed', color='black', linewidth=2.0)
    axs[1].plot( scaleVecX(-ph.hv()*omeg), scaleVecY(-ph.hv()*omeg), scaleVecZ(-ph.hv()*omeg), linestyle='solid', color='black', linewidth=4.0)
    
    
    es, phs = comptonScatter(e, ph, axs[1], plot=True)
    print"e.gamma = {}, p.hnu = {}".format(e.gamma(), ph.hv()) 
    print"es.gamma = {}, ps.hnu = {}".format(es.gamma(), phs.hv()) 
    print"Energy conservation check: incident = {}, scattered = {}".format(e.gamma()+ph.hv(),es.gamma()+phs.hv())
    #es, ps = comptonScatter(e, ph, axs[1], plot=True)
    #for sc in range(3):
    #    es, ps = comptonScatter(e, ph, axs[1], plot=True)
    plt.savefig("compton_single.png")
    #sys.exit()


    ##################################################
    # create bucket of electrons
    bucket_el = mcmc.electronBucket()
#    print "created electron bucket ({})".format( bucket_el.size() )

    #pour electrons to the bucket
    for i in range(100):
        (vx, vy, vz) = [0.0, 0.6, 0.0]#randVel(0.9) #direction on unit sphere #
        e = mcmc.electron(1.0, vx, vy, vz)
        bucket_el.push_back( e )
    print "Created electron bucket and loaded with {} electrons".format( bucket_el.size() )


    #Create photon bucket and pour photons to the bucket
    bucket_ph = mcmc.photonBucket()
    for i in range(100):
        (vx, vy, vz) = [1.0, 0.0, 0.0] #randVel(1.0) #direction on unit sphere
        EE = 1e-4 # E = h\nu/mc^2 

        ph = mcmc.photon( EE, vx, vy, vz )
        bucket_ph.push_back( ph )
    print "Created photon bucket and loaded with {} photons".format( bucket_ph.size() )
    
    
    xmin = 5.0e-5
    xmax =  5.0e-4
    axs[0].set_xlim(xmin, xmax)
    axs[0].set_ylim(0.2, 1.1)
    axs[0].minorticks_on()
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')

    
    ## constants
    deltat=1.0e5  		# timestep
    V= 1.0 * sigma_T	# considered volume
    
    for sc in range(5):
        bucket_el, bucket_ph = lp_ComptonScatter(bucket_el, bucket_ph, axs, deltat, V)
        print "Scattering # {}".format( sc+1 )
        

    plt.savefig("Compton_bucket.png")
    #plt.show()
 

