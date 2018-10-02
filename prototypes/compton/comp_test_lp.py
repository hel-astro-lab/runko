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
import el_ph

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

def randVelVec(vabs):
    (vphi, vthe) = randUnitSphere()
    vx, vy, vz = sph2cart(vabs, vphi, vthe)
     
    velVec = np.array([vx, vy, vz])
    
    return velVec


def unitVecX(vec):
    return np.array([ 0.0, vec[0] ])/np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
def unitVecY(vec):
    return np.array([ 0.0, vec[1] ])/np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
def unitVecZ(vec):
    return np.array([ 0.0, vec[2] ])/np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)

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
        jj=1.0
        sum=jj**(-3)
        while (1.202*xi1) > (sum+(jj+1)**(-3)):
            jj=jj+1.0
            sum = sum + jj**(-3)
        xi=jj+1.0
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
            zeta = - kTe * np.log( xi1 * xi2 * xi3 * xi4)
            
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
def comptonScatter(e, p, w_new, r_i_new, axs, plot):

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

#    print M    
#    print"Initial photon direction: {} {} {}".format(omega[0], omega[1], omega[2])
#
#   To add exception when beta || omega
#
#    if norm(jvec) == 0.0:
#        jvec = uCross(kvec,beta0) ## TO BE CHANGED!
#        ivec = uCross(jvec, kvec)
    cosalpha = np.dot( kvec, beta0 )    # cosine of the angle between electron and k-vector
    alpha = np.arccos(cosalpha)
    sinalpha = np.sin(alpha) # sine of the same angle
  
    ################################################## 

    mu = cosalpha  # in ijk frame angle between e and ph equals to the angle between e and k vector
#    print"mu = {}".format(mu)
        
    y = p.hv() * e.gamma() * (1.0 - mu*e.beta() ) # initial photon and electron 4-product
    
    
    #scatter
    done = False
    OmegaOmegap = 0.0
    while not(done):
        z1 = np.random.rand()
        z2 = np.random.rand()
        z3 = np.random.rand()
        
        # draw new possible angles
        mup  = (e.beta() + 2.0*z1 - 1.0)/(1.0 + e.beta()*(2.0*z1 - 1.0))  # cos(alpha') = e \dot Omega'
        phip = 2.0*np.pi*z2               # azimuthal angle calculated from (-j)
#         mup  = 1.0          ### tmp
#         phip = -np.pi/2.0   ### tmp            # azimuthal angle calculated from (-j)
 
        sinalphap = np.sqrt( 1.0 - mup**2 ) #

        OmegaOmegap = mu * mup - sinalphap *  np.sin(phip) * sinalpha  # angle between incoming and outgoing photons
#         print"OmegaOmegap = {}, mu = {}".format(OmegaOmegap,mu)
        
        yp = y / (1.0 + p.hv() * (1.0 - OmegaOmegap) / ( e.gamma() * (1.0 - mup*e.beta()) ))

        YY = yp/y + (yp/y)**3 + (yp/y)**2 *( (1.0/yp - 1.0/y)**2 - 2.0*( 1.0/yp - 1.0/y) )
 
        if YY > 2.0*z3:
            done = True
    #now we have scattered successfully

    #new energy
    hvp = yp / ( e.gamma()*(1.0 - mup*e.beta()) ) 
#     print"e. gamma = {}".format(e.gamma())
#     print"hvp/hv = {}".format(hvp/p.hv())
#     tmp = (1.0 - mu*e.beta())/(1.0 - mup*e.beta() + p.hv()*(1.0-OmegaOmegap))
#     print"direct check hvp/hv = {}".format(tmp)
    
    #hvp2 = p.hv() * (1. - e.v() * mu) / (1. - e.v() * mup + p.hv() * (1. - OmegaOmegap)/ e.gamma()) # energy test
    
    #print"hv = {}, hvp2 ={}".format(hvp,hvp2)   # compare new photon energy calculated different ways
    #print"hvp*p.hv*(1-OOp) = {}, y-yp = {}".format(hvp*p.hv()*(1.-OmegaOmegap),y-yp) # check if the quantities are conserved
    
    
    #new direction in ijk coordinate system
    Omegap_ijk = np.array( [mup*sinalpha + sinalphap * np.sin(phip) * mu, -sinalphap * np.cos(phip),
                            mup * cosalpha - sinalphap * np.sin(phip) * sinalpha] )
    
    Omegap_ijk = Omegap_ijk/norm(Omegap_ijk)   
#    print"Norm Omegap_ijk = {}".format(norm(Omegap_ijk)) 
#    print"Scattered photon direction (ijk): {} {} {}".format(Omegap_ijk[0], Omegap_ijk[1], Omegap_ijk[2])
#    print"hnup / hnu = {}".format(hvp/p.hv())
    Omegap = np.dot( np.linalg.inv(M), Omegap_ijk ) # transferring back to lab system

    phs = el_ph.photon( hvp, Omegap[0], Omegap[1], Omegap[2], r_i_new[0], r_i_new[1], r_i_new[2], w_new )

    # scattered electron parameters
    gammaes = e.gamma() + p.hv() - hvp		
    vxes = ( e.gamma() * e.vx() + p.hv() * omega[0] - hvp * Omegap[0] ) / gammaes
    vyes = ( e.gamma() * e.vy() + p.hv() * omega[1] - hvp * Omegap[1] ) / gammaes
    vzes = ( e.gamma() * e.vz() + p.hv() * omega[2] - hvp * Omegap[2] ) / gammaes
    es = el_ph.electron()
    
    es.loadVelComponents( vxes, vyes, vzes )
    ves= np.array([ vxes, vyes, vzes ])
    

    return es, phs	### ! return scattered electron and photon
#     return e, phs # return scattered photon and initial electron




## Monte-Carlo Compton scattering for bucket of electrons and photons according to
## Pozdnyakov, Sobol, Sunyaev prescription
def PSS_ComptonScatter(bucket_el, bucket_ph, axs, deltat, cell_R, tau, kTbb, kTe, xmin, xmax):
    
    nbins=100
    Nph = bucket_ph.size()
    Nel = bucket_el.size()
    e = el_ph.electron()



    hv_init = np.zeros(Nph) # initial photons energies storage
    hv_sc = np.zeros(Nph) # scattered photons storage
    hnuphnu = np.zeros(Nph) # energy ratio array
    ang = np.zeros(Nph) # testing angular distribution of scattered photons
    angel = np.zeros(Nph) # testing angular distribution of recoil electrons 
    cosalpha = np.zeros(Nph) # testing energy ratio dependence on scattering angle
    
    i=0
    w_crit = 1.e-9
    w_esc = np.zeros(Nph)
    
    for i in range(Nph):
        
        p = bucket_ph.get(i)
        if p.w() < w_crit: 
            continue
        i_it=0

        omega = np.array([ p.vx(), p.vy(), p.vz() ])         # unit vector in photon direction in lab frame
        r_i = np.array([p.x(), p.y(), p.z()])                # position of the chosen photon
        l_free = lambdaFree(p, cell_R, tau, kTe)
        l_i = - np.dot(r_i,omega) + np.sqrt( cell_R**2 - np.linalg.norm(r_i)**2 + np.dot(r_i,omega)**2 )

        P_esc = np.exp(-l_i/l_free)
        
        w_esc[i] = p.w()*P_esc     # weight of the escaping part of photons
        w_new = p.w()*(1. - P_esc) # new weight of the non-escaping photon
        
        # new scattering point
        l_ph = - l_free * np.log( 1. - np.random.rand() * (1 - P_esc) )
        r_i_new = np.array([p.x()+l_ph * omega[0], p.y()+l_ph * omega[1], p.z()+l_ph * omega[2]]) 

        #### instead of bucket, use MC sampling of Maxwell:
        omega = np.zeros(3)
        done = False
        while (not done):
            ppel = maxwellSample(kTe)      ###  Maxwellian distribution  
            eGamma = np.sqrt(ppel**2 + 1.0)
            vel = ppel/eGamma # electron velocity
            #eVel = randVelVec(vel)        ###  Homogeneous in angles
            vz = 2.*np.random.rand() - 1.0
            xi2=np.random.rand()
            vx = np.sqrt(1.0 - vz**2) * np.cos(2.*np.pi*xi2)
            vy = np.sqrt(1.0 - vz**2) * np.sin(2.*np.pi*xi2)
            vel0 = (vx, vy, vz)
        
            vrel = 1.0 - np.dot(vel0, omega) * vel
        
            xi = 2.0 * eGamma * p.hv() * vrel # product of initial electron and photon 4-momenta
            if xi < 0.5:
                s_til = 1./3. + 0.141*xi - 0.12 * xi**2 + (1. + 0.5*xi) / (1. + xi)**2
            elif xi < 3.5:
                s_til = ( np.log( 1. + xi ) +0.06 ) / xi
            else: 
                s_til = ( np.log( 1. + xi ) +0.5 - 1./(2. + 0.076 * xi) ) / xi        
            
            eta_limit = 0.375 * s_til * vrel 
            eta_rand = np.random.rand()
            if eta_rand < eta_limit: 
                done=True
        
        e.loadVelComponents(vel*vel0[0], vel*vel0[1], vel*vel0[2])               ###
        
        es, ps = comptonScatter(e, p, w_new, r_i_new, axs, True)
        omegap = np.array([ ps.vx(), ps.vy(), ps.vz() ]) 
        hv_init[i] = p.hv()*p.w()
        hv_sc[i] = ps.hv()*ps.w()+p.hv()*(ps.w()-p.w())
        bucket_ph.replace(i, ps)


    phhist, edges = np.histogram(hv*511.0, np.logspace(np.log10(xmin), np.log10(xmax), nbins), weights=w_esc)
#    axs[0].bar(edges[:-1], edges[:-1]**2*phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# number of photons per log energy
    axs[0].bar(edges[:-1], phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# number of photons per log energy
    
#     phhist, edges = np.histogram(hv*511., bins=nbins, weights=w_esc)
#     axs[0].bar(edges[:-1], 114.*phhist.astype(np.float32)/Nph, width=np.diff(edges))# number of photons per log energy

# Check of average energy change in the interactions
    ave = np.sum(hv_sc/hv)/Nph
    print"Ave nup/ nu = {}, 4.0*kTe = {}".format(ave,4.0*kTe)
#    print"Ave nup/ nu = {}, 4.0*(e.gamma()**2-1.0)/3.0 = {}".format(ave,4.0*(e.gamma()**2-1.0)/3.0)
    
#    cset = axs.contourf(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
#    cset = axs.contourf(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
#    cset = axs.contourf(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)
        
    return bucket_el, bucket_ph
    

### !!! Being tested at the moment !!!



### Calculating mean free path according to sect. 9.4 of Pozdnyakov, Sobol, Sunyaev 1983
### Not optimized: calculates only electron-dependent quantities for each photon
def lambdaFree(p, cell_R, tau, kTe):
        
    n_til = 40
    
    xi_max = 4.0 * p.hv() * (1.0 + kTe * np.log( 2.0*n_til ) ) 
    
    if xi_max < 1.0e-2:
        l_free = cell_R/tau
        return l_free
    
    g_sum=0.0
    phi_sum=0.0
    for i in range(n_til):
        gamma_til = 1.0 - kTe * np.log((i+0.5)/n_til)
        g_sum = g_sum + gamma_til * np.sqrt(gamma_til**2 -1.0)

        xi_p = 2.0 * p.hv() * (gamma_til + np.sqrt(gamma_til**2 - 1.0))
        if xi_p < 0.5:
            phi_p = xi_p**2/6.0 + 0.047*xi_p**3 - 0.03*xi_p**4 + xi_p**2 / (1. + xi_p) / 2.
        elif xi_p < 3.5:
            phi_p = (1. + xi_p) * np.log( 1. + xi_p ) - 0.94 * xi_p - 0.00925
        else: 
            phi_p = (1. + xi_p) * np.log( 1. + xi_p ) - xi_p/2.0 - 13.16 * np.log(2. + 0.076 * xi_p) + 9.214
        
        xi_m = 2.0 * p.hv() * (gamma_til - np.sqrt(gamma_til**2 - 1.0))
        if xi_m < 0.5:
            phi_m = xi_m**2/6.0 + 0.047*xi_m**3 - 0.03*xi_m**4 + xi_m**2 / (1. + xi_m) / 2.
        elif xi_m < 3.5:
            phi_m = (1. + xi_m) * np.log( 1. + xi_m ) - 0.94 * xi_m - 0.00925
        else: 
            phi_m = (1. + xi_m) * np.log( 1. + xi_m ) - xi_m/2.0 - 13.16 * np.log(2. + 0.076 * xi_m) + 9.214
        
        phi_sum = phi_sum + phi_p - phi_m

    
    g_sum = g_sum * cell_R / 0.375 /tau    
    l_free = g_sum * (2.0 * p.hv())**2 / phi_sum
    
    return l_free
    


## Monte-Carlo Compton scattering for bucket of electrons and photons
## Choose photon, compute whether the scattering occurs at all using maximal (Thomson) cross-section, 
## choose electron to interact with, compute real cross-section, do scattering (by calling comptonScatter)
def lp_ComptonScatter(bucket_el, bucket_ph, bucket_phesc,axs, deltat, cell_R, n_e, kTbb, kTe, xmin,xmax,sc, n_esc):
    
    
#    sigma_T=8.0/3.0   # Thomson cross section in units of \pi r_e^2
    
    nbins=100

#     bucket_phesc = el_ph.photonBucket() # escaping photons bucket

    Nph = bucket_ph.size()
    Nel = bucket_el.size()
    hv = np.zeros(Nph) # initial photons energies storage
    for i in xrange(Nph):
        p = bucket_ph.get(i)
        hv[i] = p.hv()
#         print"hnu_ = {}".format(p.hv())

   
    hv_sc = np.zeros(Nph) # scattered photons storage
    hnuphnu = np.zeros(Nph) # energy ratio array
    ang = np.zeros(Nph) # testing angular distribution of scattered photons
    angel = np.zeros(Nph) # testing angular distribution of recoil electrons 
    cosalpha = np.zeros(Nph) # testing energy ratio dependence on scattering angle
#    e = mcmc.electron()       # scattering electron
    sum_el_target_weights = bucket_el.size()  # sum of electron weights, now equal to number of electrons
    
#    P_it_max = sum_el_weights * 2.0 * sigma_T / cell_V  # 1/P_it_max = minimal free length (P_it_max=maximal partial interaction rate) 
    P_it_max = 4.0e3*n_e # *(sum_el_target_weights/bucket_el.size())  
                               # 1/P_it_max = minimal free length (P_it_max=maximal partial interaction rate) 
                               # P_it_max = 4e3 * (n_e/1e17) = 4e3 * sum(w*_t)
                               # sum(w*_t) - sum only over those electrons that participate in the interaction
                               
    i=0
    n_esc=0
    for i in range(Nph):
        
#         print i
        p = bucket_ph.get(i)
        hnu=p.hv()
        r_i = np.array([p.x(), p.y(), p.z()])
        rr = np.sqrt(p.x()**2 + p.y()**2 + p.z()**2)                # position of the chosen photon
        
        if rr > cell_R: ### replace the escaped photon with the seed one 
            continue
#             vz = 2.*np.random.rand() - 1.0
#             xi = np.random.rand()
#             vx = np.sqrt(1. - vz**2) * np.cos(2.*np.pi*xi)
#             vy = np.sqrt(1. - vz**2) * np.sin(2.*np.pi*xi)
#             EE = 6.4e-1/511.0 # E = h\nu/mc^2 ; Fe K alpha line
# #             EE = bbodySample(kTbb)
#             ph = el_ph.photon( EE, vx, vy, vz, 0.0, 0.0, 0.0, 1.0 )
#             bucket_ph.replace(i, ph)
#             p = bucket_ph.get(i)
#             hnu=p.hv()
#             r_i = np.array([p.x(), p.y(), p.z()])
#             rr = np.sqrt(p.x()**2 + p.y()**2 + p.z()**2)       # position of the chosen photon


        omega = np.array([ p.vx(), p.vy(), p.vz() ])         # unit vector in photon direction in lab frame

        z1 = np.random.rand()
        t_free = -np.log( z1 ) /  P_it_max  # mean free time spent between interactions
        #print"i={}, t_free={}".format(i, t_free)
        
        if deltat <= t_free :  # if timestep is less than mean free path, propagate current photon and go to the next photon
#            bucket_ph.propagate(i,t_free)  # update the position of photon
#            print"No scattering: deltat = {}, t_free = {}".format(deltat,t_free)
            hv_sc[i] = p.hv()
#             r_i_new = np.array([p.x()+t_free * omega[0], p.y()+t_free * omega[1], p.z()+t_free * omega[2]]) 
            r_i_new = np.array([p.x()+deltat * omega[0], p.y()+ deltat * omega[1], p.z()+ deltat * omega[2]]) 
            rr_new = np.sqrt(r_i_new[0]**2 + r_i_new[1]**2 + r_i_new[2]**2)
            ps = el_ph.photon( hv_sc[i], omega[0], omega[1], omega[2], r_i_new[0], r_i_new[1], r_i_new[2], 1.0 )
            bucket_ph.replace(i, ps)
            if rr_new> cell_R:        ### gather escaped photons
                bucket_phesc.push_back( ps )
                n_esc = n_esc + 1
        else:
            t_it = deltat-t_free
#             print deltat, t_free, t_it
                        
            #### Rewrite in according to the equation 7 of Stern et al 
            #### (different from the eq. below if weights vary and if there is a threshold for interaction)
            rnd_tmp = np.random.rand()
            j = int(np.floor( rnd_tmp * Nel ))  # choose the target electron 
            ####  to add: if the electron was already used, then choose another target till all targets are used (in this case propagate)
            e = bucket_el.get(j)
       
            beta0 = np.array([ e.vx(), e.vy(), e.vz() ]) / e.vmod() # unit vector in direction of electron motion in lab frame
            omega = np.array([ p.vx(), p.vy(), p.vz() ])         # unit vector in photon direction in lab frame
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
            P_it_real = s0 * (vrel / 2.0)  #  partial interaction rate relative to the maximal one
            
            
            z2 = np.random.rand()
            if z2 < P_it_real: 
#                if hnu < 1e-10:
#                    print"hnu = {}, i = {}".format(p.hv(),i)
                es, ps = comptonScatter(e, p, 1.0, r_i, axs, True)
                hnu = ps.hv()
                omegap = np.array([ ps.vx(), ps.vy(), ps.vz() ]) 
#                 bucket_el.replace(j, es) ### replace electron with the scattered one
                r_i_new = np.array([p.x()+t_free*omega[0]+t_it*omegap[0], 
                                    p.y()+t_free*omega[1]+t_it*omegap[1], 
                                    p.z()+t_free*omega[2]+t_it*omegap[2]]) 
                rr_new = np.sqrt(r_i_new[0]**2 + r_i_new[1]**2 + r_i_new[2]**2)
#                print"r_new = {}, i = {}, hnu = {}".format(r_i_new, i, hnu*511.)
                ps = el_ph.photon( hnu, omegap[0], omegap[1], omegap[2], r_i_new[0], r_i_new[1], r_i_new[2], 1.0 )

#               bucket_ph.propagate(i,t_free)  # update the position of photon for the path before interaction
                bucket_ph.replace(i, ps)
                if rr_new> cell_R:                 ### gather escaped photons
                    bucket_phesc.push_back( ps )
                    n_esc = n_esc + 1
                
                if np.abs(e.gamma()+p.hv()-es.gamma()-ps.hv()) > 1e-9:
                    print"Energy not concerved, i={}, e={}, hv={}, es={}, hvs={}, diff={}".format(i,e.gamma(),p.hv(),es.gamma(),ps.hv(),e.gamma()+p.hv()-es.gamma()-ps.hv()) 
                
                hnuphnu[i] = (1.0 - np.dot( omega, ([ e.vx(), e.vy(), e.vz() ])))/ \
                    (1.0 - np.dot( omegap, ([ e.vx(), e.vy(), e.vz() ])) + p.hv()*(1.-cosalpha[i])/e.gamma())
            else:
                r_i_new = np.array([p.x()+deltat* omega[0], p.y()+deltat* omega[1], p.z()+deltat* omega[2]]) 
                ps = el_ph.photon( p.hv(), omega[0], omega[1], omega[2], r_i_new[0], r_i_new[1], r_i_new[2], 1.0 )
                bucket_ph.replace(i, ps)
                rr_new = np.sqrt(r_i_new[0]**2 + r_i_new[1]**2 + r_i_new[2]**2)

                if rr_new > cell_R:                 ### gather escaped photons
                    bucket_phesc.push_back( ps )
                    n_esc = n_esc + 1


            
#     hv_esc = np.zeros(Nph)
#     for j in range(Nph):
#         ph = bucket_ph.get(j)
# #         print"hnu = {}".format(p.hv())
#         rr = np.sqrt(ph.x()**2 + ph.y()**2 + ph.z()**2)
# #        print rr
#         if rr > cell_R:
#             hv_esc[j]=ph.hv()
#     if n_esc > 100:
#         hv_esc = np.zeros(Nph)
#         for j in range(n_esc):
#             ph = bucket_phesc.get(j)
#             hv_esc[j]=ph.hv()
#             
#         phhist, edges = np.histogram(hv_esc*511.0, bins = 10 ** np.linspace(np.log10(xmin), np.log10(xmax), nbins))
# #         phhist, edges = np.histogram(hv_esc*511.0, np.logspace(np.log10(xmin), np.log10(xmax), nbins))
#         if np.amax(phhist) < 1e-10:
#             print"No radiation escapes"
#         elif sc > 5:
# #             axs[0].bar(edges[:-1], edges[:-1]*phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# E L_E
#             axs[0].bar(edges[:-1], phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# number of photons per log energy
#    print phhist

#    axs[0].bar(edges[:-1], edges[:-1]**2*phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# number of photons per log energy

    return bucket_el, bucket_ph, bucket_phesc, n_esc



if __name__ == "__main__":

    ## constants
    nbins = 100
    n_steps = 50
    Nel = 1000000
    Nph = Nel
    
    n_e = 1.0e0   # electron number density in units of 1e17 cm^-2
    deltat =  1.0e-4  # mean free time for Thomson scattering is 2.5e-4 * (1e17 / n_e)  (seconds)
    cell_R = 1.5e7 / 3.0e10  # radius of the sphere (in cm) divided by the speed of light
    
    
    print "tau = {}".format(2.e3*n_e*cell_R)
#     deltat=1.0e-1  		# timestep = R/c; 1 time_unit = 3.33e-4 s
#     tau=1.0e-0   # optical depth of the region (for testing, the region is a sphere)    
#     cell_V= 1.66667e-3 * (3.0/4.0/np.pi)**0.5 * (Nel * 8.0/3.0 /tau)**1.5 	# considered volume in units of (\pi r_e^2)**1.5
    
    
    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(12,6))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 2)
    
    axs = []
    axs.append( plt.subplot(gs[0, 0]) )
    plt.xlabel('hnu (keV)')
    plt.ylabel('E F_E')
#     plt.ylabel('E N_ph(E)')
    axs.append( plt.subplot(gs[0, 1]) )
    plt.xlabel('p/mc')
    plt.ylabel('p f(p)')
    
    xmin = 1.0e-2
    xmax =  2.0e+3
#     xmin = 3.0e-3
#     xmax =  1.0e+5
#     xmin = 3.5e-0
#     xmax =  11.0e+0
    axs[0].set_xlim(xmin, xmax)
    axs[0].set_ylim(1e-4, 1.2e-0)
#    axs[0].set_ylim(0., 3.0e-2)
    axs[0].minorticks_on()
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    

    elxmin = 3.e-2
    elxmax =  3.e+1
    axs[1].set_xlim(elxmin, elxmax)
    axs[1].set_ylim(1e-4, 3.0e0)
    axs[1].minorticks_on()
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')

    e = el_ph.electron()

    ##################################################
    # create bucket of electrons
    bucket_el = el_ph.electronBucket()
    e = el_ph.electron()
    kTe = 1.0e+2/511.0 # electron temperature in units of mc^2
    for i in range(Nel):
#        ppel=1.58
        ppel = maxwellSample(kTe)      ###  Maxwellian distribution  
        gammael = np.sqrt(ppel**2 + 1.0)
        vell = ppel/gammael # electron velocity
#        (vx, vy, vz) = [0.0, 1.0e-9, np.sqrt(vell**2-1.e-18)]    ### Electron beam
#        (vx, vy, vz) = [0.0, 0, vell]    ### Electron beam
#        (vx, vy, vz) = randVel(vell)        ###  monoenergetic, homogeneous in angles electron beam
        vz = 2.*np.random.rand() - 1.0
        xi = np.random.rand()
        vx = np.sqrt(1. - vz**2) * np.cos(2.*np.pi*xi)
        vy = np.sqrt(1. - vz**2) * np.sin(2.*np.pi*xi)
        e.loadVelComponents(vx*vell, vy*vell, vz*vell)               ###
        bucket_el.push_back( e )
#         print"gamma = {}, elgamma = {}".format(np.sqrt(ppel**2+1.0), e.gamma())
    print "Created electron bucket and loaded with {} electrons".format( bucket_el.size() )

    p_el = np.zeros(Nel) # initial electrons Lorentz factors storage
    for k in xrange(Nel):
        e = bucket_el.get(k)
        p_el[k] = np.sqrt(e.gamma()**2-1.0)
    
    elhist, edges = np.histogram(p_el, np.logspace(np.log10(elxmin), np.log10(elxmax), nbins))
    axs[1].bar( edges[:-1], 1.*elhist/Nel, width=np.diff(edges) )# number of electrons per log of their momenta
#     axs[7].plot( edges[:-1], 1.*elhist/Nel/np.log10(edges[1]/edges[0]) )# number of electrons per log of their momenta
#    elhist, edges = np.histogram(p_el, nbins)
    
    ### comparing the sample distribution using Maxwell-Juttner for f(p)
#    kTe=0.2
    if kTe < 1.5e-3:
        fmaxjutt =  np.sqrt(2./np.pi)*edges**2 * np.exp(-edges**2/kTe/2.)/ kTe**1.5
        axs[7].plot(edges, fmaxjutt*1.0/Nel)#, "k-")
    else:
        fmaxjutt =  edges**2 * np.exp(-np.sqrt(edges**2+1.0)/kTe) / kTe / special.kv(2, 1./kTe)
#        axs[7].plot(edges, fmaxjutt*edges*np.log(10.0))#, "k-")

    
    
    #Create photon bucket and pour photons to the bucket
    bucket_ph = el_ph.photonBucket()
    bucket_phesc = el_ph.photonBucket() # escaping photons bucket
    kTbb = 1.5e-2/511.0 ### Photon blackbody temperature in units of mc^2; 2e-4 = 0.1 keV
#    E_line=6.4e/511.0  ### keV;  Fe K aloha line
#     E_line=1e-1/511.0  ### keV;  Fe K aloha line
    for i in range(Nph):
#        (vx, vy, vz) = [1.0, 0.0, 0.0] 
        vz = 2.*np.random.rand() - 1.0
        xi = np.random.rand()
        vx = np.sqrt(1. - vz**2) * np.cos(2.*np.pi*xi)
        vy = np.sqrt(1. - vz**2) * np.sin(2.*np.pi*xi)
#         EE = E_line # E = h\nu/mc^2 
        EE = bbodySample(kTbb)
### central injection of photons
#        ph = el_ph.photon( EE, vx, vy, vz, 0.0, 0.0, 0.0, 1.0 )
### random position of photons in a unit sphere
        costh = 2.*np.random.rand() - 1.0
        phi = 2.*np.pi*xi
        r = cell_R*(2.*np.random.rand() - 1.0)
        z = r * costh
        x = r * np.sqrt(1. - costh**2) * np.cos(phi)
        y = r * np.sqrt(1. - costh**2) * np.sin(phi)
        ph = el_ph.photon( EE, vx, vy, vz, x, y, z, 1.0 )
#        print "{} {} {}".format(ph.x(), ph.y(), ph.z())
        bucket_ph.push_back( ph )
    print "Created photon bucket and loaded with {} photons".format( bucket_ph.size() )
    hv = np.zeros(Nph) # initial photons energies storage
    for i in xrange(Nph):
        p = bucket_ph.get(i)
        hv[i] = p.hv()
#         print"hnu = {}".format(p.hv())
#         print"MIN = (1-beta)/(1+beta) = {}".format((1.0-e.beta())/(1.0+e.beta()))
#         print"MAX = (1+beta)/(1-beta) = {}".format((1.0+e.beta())/(1.0-e.beta()))
    phhist, edges = np.histogram(hv*511., np.logspace(np.log10(xmin), np.log10(xmax), nbins))
#     axs[0].bar(edges[:-1], edges[:-1]*phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# E L_E
    axs[0].bar(edges[:-1], phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# number of photons per log energy
    
#    xx_nu = edges / E_line - 1.0
#    kk = 15. * np.sqrt(np.pi / 2. / kTe) / 22.
#    Inu = 1. + 1.5 * xx_nu - kk * np.absolute(xx_nu) + xx_nu**2 * (35./88./kTe - 17./22.) - kk * xx_nu * np.absolute(xx_nu)
##    Inu = 1. + (1.5 - kk) * xx_nu + (35./88./kTe - 17./22. - kk) * xx_nu**2
#    axs[0].plot(edges, Inu/100)

#     phhist, edges = np.histogram(hv*511., bins=nbins)
#     axs[0].bar(edges[:-1], 114.*phhist.astype(np.float32)/Nph, width=np.diff(edges))# number of photons per unit energy interval
    

#     print edges
#    axs[6].hist(hv, np.logspace(np.log10(xmin), np.log10(xmax), nbins), normed=True)# number of photons per log energy
#     print edges
    sum=np.sum(phhist*edges[:-1])/Nph
    print"Energy of photons: {}".format(sum)

    bbrad = 3.0e12*2.*edges**3/(np.exp(edges/kTbb)-1.0)
#    axs[6].plot(edges, edges*bbrad*np.log(10.0) )#  flux per log energy
    
    for i in xrange(3):
        np.random.rand()
    
    n_esc_tot=0
    n_esc=0
    for sc in range(n_steps):
#         bucket_el, bucket_ph = PSS_ComptonScatter(bucket_el, bucket_ph, axs, deltat, cell_R, tau, kTbb, kTe, xmin,xmax)
#         print "Scattering # {}".format( sc+1 )
        bucket_el, bucket_ph, bucket_phesc, n_esc = lp_ComptonScatter(bucket_el, bucket_ph, bucket_phesc, axs, deltat, cell_R, n_e, kTbb, kTe, xmin,xmax,sc, n_esc)
        print "Time step # {}".format( sc+1 )
        print "{} photons escaped".format(n_esc)  
        n_esc_tot = n_esc_tot + n_esc  
        if n_esc_tot == bucket_ph.size():
            break

    print"Sizeof bucket = {}, n_esc_tot = {}".format(bucket_ph.size(), n_esc_tot)
    hv_esc = np.zeros(n_esc_tot)
    for j in range(n_esc_tot):
        ph = bucket_phesc.get(j)
        hv_esc[j]=ph.hv()
            
#     phhist, edges = np.histogram(hv_esc*511.0, bins = np.linspace(xmin, xmax, nbins))
    phhist, edges = np.histogram(hv_esc*511.0, bins = 10 ** np.linspace(np.log10(xmin), np.log10(xmax), nbins))
#     phhist, edges = np.histogram(hv_esc*511.0, np.logspace(np.log10(xmin), np.log10(xmax), nbins))
    axs[0].bar(edges[:-1], 0.5*(edges[:-1]+edges[1:])*phhist.astype(np.float32)/n_esc_tot, width=np.diff(edges), log=True)# E L_E
#     axs[0].bar(edges[:-1], phhist.astype(np.float32)/n_esc_tot, width=np.diff(edges))# number of photons per log energy
#     axs[0].bar(edges[:-1], phhist.astype(np.float32)/n_esc_tot, width=np.diff(edges), log=True)# number of photons per log energy
#    print phhist
    axs[0].plot(edges[:-1], 3.0e-3*edges[:-1]**(0.1))
#         for i in xrange(Nph):
#             p = bucket_ph.get(i)
#             hv[i] = p.hv()
#             e = bucket_el.get(i)
#             p_el[i] = np.sqrt(e.gamma()**2-1.0)
# #             print"i={}; p_el,sc={}".format(i,p_el[i])
# 
#         elhist, edges = np.histogram(p_el, np.logspace(np.log10(elxmin), np.log10(elxmax), nbins))
# #         axs[7].plot( edges[:-1], 1.*elhist/Nel/np.log10(edges[1]/edges[0]) )# number of electrons per log of their momenta
# #        elhist, edges = np.histogram(p_el, nbins)
#         axs[1].bar( edges[:-1], 1.*elhist/Nel, width=np.diff(edges) )#, "k-")

#        phhist, edges = np.histogram(hv, np.logspace(np.log10(xmin), np.log10(xmax), nbins))
#         axs[6].plot(edges[:-1], 1.5e13*edges[:-1]*phhist*(kTbb**3)/np.log10(edges[1]/edges[0])/Nph )#, "k-")
#         axs[6].plot(edges[:-1], 1.5e10*phhist*(kTbb**3)/np.log10(edges[1]/edges[0])/Nph )#, "k-")
#        axs[0].bar(edges[:-1], phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# number of photons per log energy
#        print 1.5e13*edges[:-1]*phhist*(kTbb**3)/(np.log(edges[1]/edges[0]))/Nph
#        print elhist/(1.*Nel)

    plt.savefig("Compton_lp.png")
    #plt.show()
 

