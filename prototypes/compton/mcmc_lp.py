import numpy as np
import os, sys
import copy
#from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
#import palettable as pal
from matplotlib import cm
from astroML.plotting import hist

from scipy import stats, special
#from visualize import *
import mcmc

#set seed so we get predictable results
np.random.seed( 12 )
from mpl_toolkits.mplot3d import Axes3D


#Profiling / Timing
sys.path.insert(0, '../tools')
from timer import Timer

#sigma_T = 8.0/3.0    # Thomson cross-section = 8 \pi r_e^2 / 3, in units of \pi r_e^2


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

#######################################

#### Sample from blackbody distribution
def bbodySample(kTbb):
    
    xi1=np.random.rand()
    xi2=np.random.rand()
    xi3=np.random.rand()
    xi4=np.random.rand()
    
    #print xi1
    if 1.202*xi1 < 1: 
        xi=1.0
        #print xi
    else:
        jj=2.0
        sum=1.0 + jj**(-3)
        while (1.202*xi1 <= sum) & (1.202*xi1 > sum+(jj+1)**(-3)):
            jj=jj+1.0
            sum = sum + jj**(-3)
        xi=jj
        #print xi
        
    hv = - kTbb * np.log(xi2*xi3*xi4) / xi
  
    return hv


### Sample from Maxwellian distribution
def maxwellSample(kTe):

    done = False
    if kTe < 0.29:
        while (not done):
            xi1=np.random.rand()
            xi2=np.random.rand()
        
            xip = -1.5 * np.log(xi1)
        
            xi_limit = 0.151 * (1. + kTe * xip)**2 * xip * (2. + kTe * xip) * xi1
            if xi2**2 < xi_limit:
                pel = np.sqrt( kTe * xip * (2. + kTe * xip) )
                done = True

    else:
        while (not done):
            xi1=np.random.rand()
            xi2=np.random.rand()
            xi3=np.random.rand()
            xi4=np.random.rand()
            
            eta = - kTe * np.log( xi1 * xi2 * xi3)
            eta = - kTe * np.log( xi1 * xi2 * xi3 * xi4)
            
            if (zeta**2 - eta**2) > 1.0:
                pel = eta
                done = True

    return pel                 


########################################

# Monte Carlo Compton scattering according to Sobol 1977
# The reference frame where the scattering is calculated is chosen as:
# k = Omega (direction of incoming photon)
# j = k x beta (incoming electron is in the (i,k) plane)
# i - complementing the right-handed basis
def comptonScatter(e, p, axs, plot):

    # notes about units
    # p.hv = hv/m_e c^2
    # e.v = v/c 

    beta0 = np.array([ e.vx(), e.vy(), e.vz() ]) / e.vmod() # unit vector in direction of electron motion in lab frame
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
#        jvec = uCross(kvec,beta0) ## TO BE CHANGED!
#        ivec = uCross(jvec, kvec)
    cosalpha = np.dot( kvec, beta0 )    # cosine of the angle between electron and k-vector
    sinalpha = np.sqrt( 1.0 - cosalpha**2 ) # sine of the same angle
  
    ################################################## 

    mu = cosalpha  # in ijk frame angle between e and ph equals to the angle between e and k vector
        
    y = p.hv() * e.gamma() * (1.0 - mu*e.vmod() ) # initial photon and electron 4-product
    
    #scatter
    done = False
    OmegaOmegap = 0.0
    while not(done):
        z1 = np.random.rand()
        z2 = np.random.rand()
        z3 = np.random.rand()
        
        # draw new possible angles
        mup  = (e.vmod() + 2.0*z1 - 1.0)/(1.0 + e.vmod()*(2.0*z1 - 1.0))  # cos(alpha') = k \dot Omega'
        phip = 2.0*np.pi*z2               # azimuthal angle calculated from (-j)
        sinalphap = np.sqrt( 1.0 - mup**2 ) #

        OmegaOmegap = mu * mup - sinalphap *  np.sin(phip) * sinalpha  # angle between incoming and outgoing photons
        
        yp = y / (1.0 + p.hv() * (1.0 - OmegaOmegap) / ( e.gamma() * (1.0 - mup*e.vmod()) ))

        YY = yp/y + (yp/y)**3 + (yp/y)**2 *( (1.0/yp - 1.0/y)**2 - 2.0*( 1.0/yp - 1.0/y) )
 
        if YY > 2.0*z3:
            done = True
    #now we have scattered successfully

    #new energy
    hvp = yp / ( e.gamma()*(1.0 - mup*e.vmod()) ) 
    
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
    es = mcmc.electron()
    
    es.loadVelComponents( vxes, vyes, vzes )
    ves= np.array([ vxes, vyes, vzes ])
    
#### Plotting 3D vectors of scattered particles
    if plot:
        axs[3].plot( scaleVecX(hvp*Omegap), scaleVecY(hvp*Omegap), scaleVecZ(hvp*Omegap), alpha=0.2, linestyle='solid', color='red')#, label='Omegap' )
        axs[3].plot( scaleVecX(ves), scaleVecY(ves), scaleVecZ(ves), alpha=0.2, linestyle='dashed', color='blue')#, label='Omegap' )

    return es, phs	### ! return scattered electron and photon
#    return e, phs # return scattered photon and initial electron



# Lorentz transformation to electron rest frame (sign=1) or back to lab frame (sign=-1)
# Lorentz transformation of 4-vector vec1 to the electron frame with 4-vector vec2
# def lorentz(vec1, vec2, intt)
#    
# 	gamma = vec2[0]
# 	eph =  vec1[0]  
#    
# 	vec1spacial = (vec1[1], vec1[2], vec1[3])
# 	vec2spacial = (vec2[1], vec2[2], vec2[3])
#    
# 	pv = matmul(np.transpose(vec1spacial),vec2spacial)
# 	t = (gamma - 1.0) * pv - eph * np.sqrt(gamma**2 - 1.0)
#    
# 	eph = gamma*eph + pv * np.sqrt(gamma**2 - 1.0)
# 	vec1spacial = vec1spacial + t * vec2spacial
#    
# 	vec1 = (eph, vec1spacial[0], vec1spacial[1], vec1spacial[2])
# 
# 	return vec1
    
    
    
    
#     
# # Transformation (rotation) of vec2 in the reference frame connected with vec1 to lab frame
# def transfToLab(vec1, vec2)
#     
#     t = np.sqrt( ve1[0]**2 + vec1[1]**2 )
#         
#     if t == 0.0:
#         mtransf = np.array([0., 0., -1.], [0., 1., 0.], [1., 0., 0.]) 
#     else:
#         a=vec[0]
#         b=vec[0]
#         c=vec[0]
#         mtransf = np.array(a, -a*c/t, b/t], [b, -b*c/p, -a/p], [c, p, 0.]) 
#         
#     vec2 = np.matmul(mtransf,vec2)
#                     
#     return vec2





# Monte Carlo Compton scattering by Boris
def comptonScatterBoris(e, p, ax, plot):
        
    lorentz(p,e,1) # transformation of the photon vector from lab to electron rest frame
    
    # perform Compton scattering in the electron rest frame
    #    

    omega0 # unit vector in the direction of the incoming photon
    omegap # unit vector in the direction of the outgoing photon
    vel= np.array( vxes, vyes, vzes )
    
    omegap=transfToLab(omega0,omegap) # transformation (rotation) of the scattered photon vector to lab frame
    vel=transfToLab(omega0,vel) # transformation (rotation) of the scattered electron vector to the lab frame
    
    es.loadVel(vel)
    phs = mcmc.photon( hvp, omegap[0], omegap[1], omegap[2] )

    phs=lorentz(phs,e,-1) # transformation of the photon vector from the incoming electron rest frame to lab frame
    es=lorentz(es,e,-1) # transformation of the scattered electron vector from incoming electron rest frame to lab frame
    

    return es, phs	### ! return scattered electron and photon
#   return e, phs # return scattered photon and initial electron



### !!! Being tested at the moment !!!

## Monte-Carlo Compton scattering for bucket of electrons and photons
## Choose photon, compute whether the scattering occurs at all using maximal (Thomson) cross-section, 
## choose electron to interact with, compute real cross-section, do scattering (by calling comptonScatter)
def lp_ComptonScatter(bucket_el, bucket_ph, axs, deltat, cell_V):
    
    sigma_T=8.0/3.0   # Thomson cross section in units of \pi r_e^2
    
    Nph = bucket_ph.size()
    Nel = bucket_el.size()
    hv = np.zeros(Nph) # initial photons energies storage
    for i in xrange(Nph):
        p = bucket_ph.get(i)
        hv[i] = p.hv()
        
    xmin = 5.0e-6
    xmax =  5.0e-1
    
    phhist, edges = np.histogram(hv, np.logspace(np.log10(xmin), np.log10(xmax), 50))
    #hist = 1.0 * hist / hist.max()
    #print"hist.max = {}".format(hist.max())
    axs[6].plot(edges[:-1], edges[:-1]*phhist )#, "k-")

    kTbb=2e-4
    bbrad = 5.0e13*2.*edges**3/(np.exp(edges/kTbb)-1.0)
    axs[6].plot(edges, edges*bbrad )#, "k-")

    elgam = np.zeros(Nph) # initial electrons Lorentz factors storage
    for k in xrange(Nel):
        e = bucket_el.get(k)
        elgam[k] = e.gamma()
#         print"elgam = {}".format(elgam[k])
    
    elxmin = 1.0e-3
    elxmax =  5.0e+0

    elhist, edges = np.histogram(elgam, np.logspace(np.log10(elxmin), np.log10(elxmax), 50))
    axs[7].plot(edges[:-1], elhist )#, "k-")
        
    ### checking the sample distribution using Maxwell-Juttner formula
    kTe=0.2
    fmaxjutt = edges**2 * np.sqrt(1.0-edges**(-2)) * np.exp(-edges/kTe) / kTe / special.kv(2, 1./kTe)
    axs[7].plot(edges-1.0, fmaxjutt )#, "k-")
#    print elgam
   
    hv_sc = np.zeros(Nph) # scattered photons storage
    hnuphnu = np.zeros(Nph) # energy ratio array
    ang = np.zeros(Nph) # testing angular distribution of scattered photons
    angel = np.zeros(Nph) # testing angular distribution of recoil electrons 
    cosalpha = np.zeros(Nph) # testing energy ratio dependence on scattering angle
#    e = mcmc.electron()       # scattering electron
    sum_el_weights = bucket_el.size()  # sum of electron weights, now equal to number of electrons
    
    P_it_max = 1e-17*sum_el_weights * 2.0 * sigma_T / cell_V  # 1/P_it_max = minimal free length (P_it_max=maximal partial interaction rate) 
    ### factor 1e17 is for testing purposes only
    
    print"P_it_max = {}".format(P_it_max)
    
                                                        # cell_V-space volume; sigma_T - Thomson cross-section 
    #print"P_max = {}, bucket size = {}".format(P_it_max,bucket_ph.size())
    
    
    i_it=0    
    i=0
    
    while i < Nph:
        
        p = bucket_ph.get(i)
#        hv_sc[i] = p.hv()

        z1 = np.random.rand()
        t_free = -np.log( z1 ) /  P_it_max  # mean free time spent between interactions

        
        if deltat <= t_free :  # if timestep is less than mean free path, propagate current photon and go to the next photon
#            bucket_ph.propagate(i,t_free)  # update the position of photon
#            print"No scattering: deltat = {}, t_free = {}".format(deltat,t_free)
            hv_sc[i] = p.hv()

#            ang[i]=np.arctan2(p.vz(),np.sign(ps.vx())*np.sqrt(p.vx()**2+p.vy()**2)) # to test angular distribution
#            ang[i]=np.arctan2(p.vz(),ps.vx()) # to test angular distribution

#            print"ang [{}] = {}".format(i,ang[i])
            i=i+1
            continue
        else:
            t_it = deltat-t_free
#            print"Scattering: deltat = {}, t_free = {}".format(deltat,t_free)

            #### Rewrite in according to the equation 7 of Stern et al (different from the eq. below if weights vary)
            j = int(np.floor( np.random.rand() * Nel ))  # choose the target electron 
            
            j = i  ### for testing only: each electron interacts witch each photon
            #print"j = {}".format(j)
            ####  ??? What happens if the electron was already "used"???
            e = bucket_el.get(j)
       
            beta0 = np.array([ e.vx(), e.vy(), e.vz() ]) / e.vmod() # unit vector in direction of electron motion in lab frame
            omega = np.array([ p.vx(), p.vy(), p.vz() ])         # unit vector in photon direction in lab frame
#           mu = np.dot(beta0, omega)    # cosine of angle between incident photon and electron
            mu = np.dot(beta0, omega)    # cosine of angle between incident photon and electron

           #calculate real cross-section and real probability of interaction

            xi = e.gamma() * p.hv() * (1.0 - e.vmod() * mu) # product of initial electron and photon 4-momenta
           
            if xi < 0.01: # Approximate Taylor expansion for Compton total cross-section if xi<0.01, error about 5e-6
                s0 = 1.0 - 2.0 * xi + 5.2 * xi**2 - 9.1 * xi**3 
                # + 1144.0 * xi**4 / 35.0  - 544.0 * xi**5 / 7. + 1892.0 * xi**6 / 21.0   # higher-order terms, not needed if xi<0.01 
            else:    # Exact formula for Klein - Nishina cross section in units of sigma_T
                s0 = 3.0 / 8.0 / xi**2 * ( 4.0 + (xi - 2.0 - 2.0 / xi) * np.log(1.0 + 2.0 * xi) + 
                                       2.0 * xi**2 * (1.0 + xi) / (1.0 + 2.0 * xi)**2 )
        
            vrel = 1.0 - e.vmod() * mu  # relative velocity ("v_rel" in Stern et al. 1995), electron velocity in units of c
            P_it_real = s0 * vrel / 2.0 #  partial interaction rate relative to the maximal one
#            print"P_it_real = {}".format(P_it_real)
            
            z2 = np.random.rand()
            if z2 < P_it_real: 
#                axs[1].plot( unitVecX(-beta0), unitVecY(-beta0), unitVecZ(-beta0), linestyle='dashed', color='black', linewidth=2.0)
#                axs[1].plot( scaleVecX(-ph.hv()*omeg), scaleVecY(-ph.hv()*omeg), scaleVecZ(-ph.hv()*omeg), linestyle='solid', color='black', linewidth=4.0)

                es, ps = comptonScatter(e, p, axs, True)
                omegap = np.array([ ps.vx(), ps.vy(), ps.vz() ]) 
                cosalpha[i] = np.dot( omega, omegap )
#!!! Commented out for testing purposes; should be returned!  #            bucket_el.replace(j, es)
                angel[i]=np.arctan2(es.vx(),es.vz()) # to test recoil electrons angular distribution 

#               bucket_ph.propagate(i,t_free)  # update the position of photon for the path before interaction
                bucket_ph.replace(i, ps)
                hv_sc[i] = ps.hv()
                if np.abs(e.gamma()+p.hv()-es.gamma()-ps.hv()) > 1e-10:
                    print"Ener,gy not concerved, i={}, e={}, hv={}".format(i,e.gamma(),p.hv()) 
#                print"e.gamma = {}, p.hnu = {}".format(e.gamma(), p.hv()) 
#                print"es.gamma = {}, ps.hnu = {}".format(es.gamma(), ps.hv()) 
#                print"Energy conservation check: incident = {}, scattered = {}".format(e.gamma()+p.hv(),es.gamma()+ps.hv())
#               bucket_ph.propagate(i,t_it)  # update the position of the scattered photon
#               bucket_el.propagate(j,t_it)  # update the position of the scattered photon
               #print"# {} scattered successfully, new energy = {}".format(i, hv_sc[i])

#                ang[i]=np.arctan2(ps.vz(),np.sign(ps.vx())*np.sqrt(ps.vx()**2+ps.vy()**2)) # to test angular distribution 
                ang[i]=np.arctan2(ps.vz(),ps.vx()) # to test photon angular distribution 
                
                hnuphnu[i] = (1.0 - np.dot( omega, ([ e.vx(), e.vy(), e.vz() ])))/ \
                    (1.0 - np.dot( omegap, ([ e.vx(), e.vy(), e.vz() ])) + hv[i]*(1.-cosalpha[i])/e.gamma())


#                print"ang [{}] = {}".format(i,ang[i])
                i_it=i_it+1
                i=i+1
#                print"Scattering succeeded, iterations = {}".format(i_it)
                i_it=0
            else:
                i_it=i_it+1 
#                print"No scattering, i = {}, i_it = {}".format(i, i_it)
#                i=i-1      # if the scattering did not succeed with the chosen electron,
                continue   # try simulating scattering again with the same photon
                
    nbins=100
     
### Photon energy distribution
    #print"hist = {}".format(hist)
    #hist = 1.0 * hist / hist.max()
    ehist, edges = np.histogram(hv_sc, np.logspace(np.log10(xmin), np.log10(xmax), nbins), normed=True) #, density=True)
    axs[6].plot(edges[:-1], edges[:-1]*ehist )
    #hist(hv_sc, bins='blocks', log=True, normed=True)
    #print hist
    #axs[0].plot(edges[:-1], ehist )
 
    theta=np.zeros(nbins)
    anghist, ang_edges = np.histogram(ang, bins=nbins, range=(-np.pi, np.pi) )
    angelhist, angel_edges = np.histogram(angel, bins=nbins, range=(-np.pi, np.pi) )
#    print"anghist = {}".format(anghist)
#    print"ang_edges = {}".format(ang_edges)
    for jj in xrange(nbins):
        theta[jj]= 0.5*(ang_edges[jj]+ang_edges[jj+1])
#    print"theta = {}".format(theta)
    axs[1].plot( theta, anghist )    # distribution of scattered photons in polar coordinates
    axs[4].plot( theta, angelhist )  # distribution of recoil electrons in polar coordinates

### Scattered photons distribution histogram
#    axs[2].plot(theta, 3.0*anghist/(anghist[nbins/2-1]+anghist[nbins/2]+anghist[nbins/2+1]) )
    axs[2].plot(theta, anghist )
### analytical function for differential cross-section d \sigma (\alpha) / d \Omega    
    sigma_compare = 0.5*(1.0 + np.cos(theta)**2) * (1.0 + hv[0]*(1.0 - np.cos(theta)))**(-2) \
                    * (1.0 + hv[0]**2 * (1.0 - np.cos(theta))*2 / (1.0 + np.cos(theta)**2) \
                    /(1.0 + hv[0]*(1.0 - np.cos(theta))))
    axs[2].plot(theta, sigma_compare ) 
    
### Recoil electrons distribution    
    axs[5].plot(theta, 3.0*angelhist/(angelhist[nbins/2-1]+angelhist[nbins/2]+angelhist[nbins/2+1]) )

### Distribution of energy increase, \nu'/\nu, with respect to the scattering angle
    #print np.ndarray.shape(cosalpha)
    #print np.ndarray.shape(hv_sc/hv)
    #sys.exit()
    #cosalphahist, cosalphaedges = np.histogram(cosalpha, bins=nbins, range=[-1., 1.], \
    #                                           weights=hv_sc/hv, density=True)
    #axs[6].plot(cosalphaedges[:-1], cosalphahist)
    statistic, bin_edges, binnumber =  stats.binned_statistic(cosalpha,hv_sc/hv, statistic='mean', bins=nbins)
    axs[0].plot(np.arccos(bin_edges[:-1]), statistic)

    #hnuphnu = 1./(1.+ hv*(1.-cosalpha))
    statistic, bin_edges, binnumber =  stats.binned_statistic(cosalpha, hnuphnu, statistic='mean', bins=nbins/2)
    axs[0].plot(np.arccos(bin_edges[:-1]), statistic)



# Check of average energy change in the interactions
    ave = np.sum(hv_sc/hv)/Nph - 1.0
    print"Ave Delta nu/ nu = {}, rhs = {}".format(ave,4.0*(e.gamma()**2-1.0)/3.0)
    
#    cset = axs.contourf(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
#    cset = axs.contourf(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
#    cset = axs.contourf(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)
        
    return bucket_el, bucket_ph
    
    

# def visualize(ax, slab, params):
# 
#     xs = slab.xloc
#     ys = slab.yloc
#     zs = slab.zloc
# 
#     ax.plot(xs, ys, zs, 
#             color='black',
#             marker='.'
#             )
#             
# 
#     return 



if __name__ == "__main__":

    ## constants
    Nel=100
    Nph=Nel
    deltat=1.0e0  		# timestep
    cell_V= 1.0e0*(8.0/3.0)**1.5 	# considered volume in units of (\pi r_e^2)**1.5
    
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')

    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(9,9))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(3, 3)
    gs.update(hspace = 0.5)
    plt.subplots_adjust(wspace=0.38)
    
    axs = []
    axs.append( plt.subplot(gs[0, 0]) )
    plt.xlabel('alpha')
    plt.ylabel('E_sc/E_init (alpha) ')
    axs.append( plt.subplot(gs[0, 1], projection='polar') )
    axs.append( plt.subplot(gs[0, 2]) )
    plt.xlabel('alpha from initial photon direction')
    plt.ylabel('sigma(alpha) / sigma(0)')
    axs.append( plt.subplot(gs[1, 0], projection='3d') )
    plt.xlabel('x')
    plt.ylabel('y')
    axs.append( plt.subplot(gs[1, 1], projection='polar') )
    axs.append( plt.subplot(gs[1, 2]) )
    plt.xlabel('theta')
    plt.ylabel('recoil electrons')
    axs.append( plt.subplot(gs[2, 0]) )
    plt.xlabel('hnu/mc^2')
    plt.ylabel('E L_E')
    axs.append( plt.subplot(gs[2, 1]) )
    plt.xlabel('gamma/mc^2')
    plt.ylabel('f_gamma')
    
#    plt.zlabel('z')
#    axs.append( plt.subplot(gs[2]) )
    

    ################################################## 
    # Normalize 3d plot to make aspect ratio unity
    MAX = 1.0
    for direction in (-1, 1):
        for point in np.diag(direction * MAX * np.array([1,1,1])):
            #axs[1].plot([point[0]], [point[1]], [point[2]], 'w')
            #axs[2].plot([point[0]], [point[1]], [point[2]], 'w')

            axs[6].plot([point[0]], [point[1]], [point[2]], 'w')
#            axs[2].plot([point[0]], [point[1]], [point[2]], 'w')

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
#     print "Compton scatter single interaction..."

    e = mcmc.electron()
    verand=0.8
    gamme=5.0005
    vzz=np.sqrt(1.-gamme**(-2))
    (vx, vy, vz) = [0.0, 0.0, vzz] #randVel(verand) #(0.1, 0.0, 0.0) #randVel(0.5) 
    e.loadVelComponents(vx, vy, vz)
#     print "target electron with beta: {} gamma: {}".format(e.beta(), e.gamma() )

    (vx, vy, vz) = [1.0, 0.0, 0.0] #randVel(1.0) #direction on unit sphere (-1.0, 0.0, 0.0) 
    EE = 1e-3  # in units of mc2
    ph = mcmc.photon( EE, vx, vy, vz )

    #visualize starting point of collision
    beta0 = np.array([ e.vx(), e.vy(), e.vz() ]) / e.vmod()  # unit vector in electron direction
    
    omeg = np.array([ ph.vx(), ph.vy(), ph.vz() ])         #photon unit vector
#    axs[1].plot( unitVecX(-beta0), unitVecY(-beta0), unitVecZ(-beta0), linestyle='dashed', color='black', linewidth=2.0)
    axs[3].plot( unitVecX(-beta0)*e.vx(), unitVecY(-beta0)*e.vy(), unitVecZ(-beta0)*e.vz(), linestyle='dashed', color='black', linewidth=2.0)
    axs[3].plot( scaleVecX(-ph.hv()*omeg), scaleVecY(-ph.hv()*omeg), scaleVecZ(-ph.hv()*omeg), linestyle='solid', color='black', linewidth=4.0)
    
    
    es, phs = comptonScatter(e, ph, axs, plot=True)
    print"e.gamma = {}, p.hnu = {}".format(e.gamma(), ph.hv()) 
    print"e.gamma^2-1 = {}".format(e.gamma()**2-1.0) 
    print"beta = {}".format(e.beta()) 
    print"es.gamma = {}, ps.hnu = {}".format(es.gamma(), phs.hv()) 
    print"Energy conservation check: incident = {}, scattered = {}".format(e.gamma()+ph.hv(),es.gamma()+phs.hv())
   #es, ps = comptonScatter(e, ph, axs[1], plot=True)
    #for sc in range(3):
    #    es, ps = comptonScatter(e, ph, axs[1], plot=True)
    #plt.savefig("compton_single.png")

    ##################################################
    # create bucket of electrons
    bucket_el = mcmc.electronBucket()
#    print "created electron bucket ({})".format( bucket_el.size() )
    e = mcmc.electron()
    kTe = 0.2 # electron temperature in units of mc^2
    #pour electrons to the bucket
    for i in range(Nel):
#         ppel=1.0653
#         gammael = np.sqrt(ppel**2 + 1.0)
#         vell = ppel/gammael
        #(vx, vy, vz) = [0.0, 0.0, vell]    ### Electron beam
        #(vx, vy, vz) = randVel(vell)       ### direction on unit sphere #
        #e.loadVelComponents(vx, vy, vz)   ### 
        
        ppel = maxwellSample(kTe)                                ###   
        gammael = np.sqrt(ppel**2 + 1.0)
        vell = ppel/gammael
        (vx, vy, vz) = randVel(vell)                    ###  Maxwellian distribution
        e.loadVelComponents(vx, vy, vz)               ###
        bucket_el.push_back( e )
#         print"gamma = {}, elgamma = {}".format(np.sqrt(ppel**2+1.0), e.gamma())
        

    print "Created electron bucket and loaded with {} electrons".format( bucket_el.size() )
    
#     eltest = mcmc.electron()
#     for k in range(Nel):
#         eltest = bucket_el.get(k)
#         elgamma = eltest.gamma()
#         print"elgamma = {}, pel = {}".format(elgamma, ppel)
#     sys.exit()
    
    #Create photon bucket and pour photons to the bucket
    bucket_ph = mcmc.photonBucket()
    kTbb = 2.e-4 ### Photon blackbody temperature in units of mc^2
    for i in range(Nph):
        #(vx, vy, vz) = [1.0, 0.0, 0.0] 
        (vx, vy, vz) = randVel(1.0) #direction on unit sphere
        #EE = 1e-6 # E = h\nu/mc^2 
        EE_new=bbodySample(kTbb) #*np.exp(-np.random.rand())
        ph = mcmc.photon( EE_new, vx, vy, vz )
        bucket_ph.push_back( ph )
    print "Created photon bucket and loaded with {} photons".format( bucket_ph.size() )
    
    
    #sys.exit()
    
    xmin = 5.0e-6
    xmax =  5.0e-1
    axs[6].set_xlim(xmin, xmax)
    axs[6].set_ylim(1e-4, 2.1e0)
    axs[6].minorticks_on()
    axs[6].set_xscale('log')
    axs[6].set_yscale('log')

    elxmin = 1.0e-3
    elxmax =  5.0e+0
    axs[7].set_xlim(elxmin, elxmax)
    axs[7].set_ylim(1e-2, 2.1e1)
    axs[7].minorticks_on()
    axs[7].set_xscale('log')
    axs[7].set_yscale('log')

#     axs[2].set_xlim(xmin, xmax)
#     axs[2].set_ylim(1e-2, 1.1)
#     axs[2].minorticks_on()
#     axs[2].set_xscale('log')
#     axs[2].set_yscale('log')


    for i in xrange(3):
        np.random.rand()
    
    for sc in range(1):
        bucket_el, bucket_ph = lp_ComptonScatter(bucket_el, bucket_ph, axs, deltat, cell_V)
        print "Scattering # {}".format( sc+1 )
        

    plt.savefig("Compton_bucket.png")
    #plt.show()
 

