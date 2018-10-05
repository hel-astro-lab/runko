import numpy as np




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


def randab(a, b):
    return a + (b-a)*np.random.rand()

# Cartesian 3D location
def rand3Dloc(conf):
    x = randab(conf.xmin, conf.xmax)
    y = randab(conf.ymin, conf.ymax)
    z = randab(conf.zmin, conf.zmax)

    return [x,y,z]


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

def rand3Dvel(vabs):
    (vphi, vthe) = randUnitSphere()
    vx, vy, vz = sph2cart(vabs, vphi, vthe)
    return [vx, vy, vz]






