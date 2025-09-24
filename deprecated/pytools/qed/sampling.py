# -*- coding: utf-8 -*- 

import numpy as np



# Sample from blackbody distribution
def sample_blackbody(kTbb):
    
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


# Sample from Maxwellian distribution
def sample_maxwellian(kTe):
    pel = 0.0
    if kTe < 0.29:
        while True:
            xi1=np.random.rand()
            xi2=np.random.rand()
        
            xip = -1.5 * np.log(xi1)
        
            xi_limit = 0.151*(1. + kTe*xip)**2*xip*(2. + kTe*xip)*xi1
            if xi2**2 < xi_limit:
                pel = np.sqrt( kTe*xip*(2. + kTe*xip) )
                break
    else:
        while True:
            xi1=np.random.rand()
            xi2=np.random.rand()
            xi3=np.random.rand()
            xi4=np.random.rand()
            
            eta  = -kTe*np.log( xi1*xi2*xi3)
            zeta = -kTe*np.log( xi1*xi2*xi3*xi4)
            
            if (zeta**2 - eta**2) > 1.0:
                pel = eta
                break
    return pel



def randab(a, b):
    return a + (b-a)*np.random.rand()

# Cartesian 3D location
def rand_3D_loc(conf):
    x = randab(conf.xmin, conf.xmax)
    y = randab(conf.ymin, conf.ymax)
    z = randab(conf.zmin, conf.zmax)

    return [x,y,z]


#random 3D direction in terms of angles phi and theta
def rand_unit_sphere():
    vphi = randab( 0.0, 2.0*np.pi ) #azimuth
    vthe = randab(-np.pi/2.0, np.pi/2.0 ) #latitude
    return (vphi, vthe)

def sph2cart(vr, vphi, vthe):
    vx = vr * np.sin(vthe) * np.cos(vphi)
    vy = vr * np.sin(vthe) * np.sin(vphi)
    vz = vr * np.cos(vthe)
    return vx, vy, vz

def rand_3D_vel(vabs):
    (vphi, vthe) = rand_unit_sphere()
    vx, vy, vz = sph2cart(vabs, vphi, vthe)
    return [vx, vy, vz]



