import numpy as np

dx = 1.0
dt = 0.1
nx = 200
nv = 100
ntime = 2000

#species
ns = 2
wp   = np.array([ 0.707107,  0.707107]) #plasma frequency
qm   = np.array([-1.000000, -1.000000]) #charge-to-mass ratio
vt   = np.array([ 1.000000,  1.000000]) #thermal velocity
vd   = np.array([ 2.000000,  -2.000000]) #drift velocity

vmax = np.array([ 9.999990,  9.999990])
vmin = np.array([-5.000000, -5.000000])


dv = np.zeros(ns)
qn = np.zeros(ns)


#spatial array sections
nxHalo = 3 #width of halo region
nxfull = nx + 2*nxHalo

xLb   = np.arange(0, nxHalo) #left boundary halo region
xLe   = np.arange(nxHalo, nxHalo + nxHalo) #left edge of computational regime

xRb   = np.arange(-nxHalo, 0) #right boundary halo region
xRe   = np.arange(-2*nxHalo, -nxHalo) #right edge of computational regime

xmid  = np.arange(nxHalo, nxfull - nxHalo) #mid computational regime
xfull = np.arange(0, nxfull) #full x array range


#velocity array sections
nvHalo = 3
nvfull = nv + 2*nvHalo
vfull = np.arange(0, nvfull) 



#--------------------------------------------------
#noise  3 - white noise
#noise = 3 


nmode = 1.0 #number of noise modes
pamp = 0.01 #positive noise amplitude
namp = 0.00 #negative noise amplitude
#pphs = 0.00
#nphs = 0.00


#diagnostics
#nplot    = 1000.000000
#emax     = 1.000000
#diagtype = [10.000000, 11.000000, 2.000000, 9.000000, 5.000000, 1.000000 ]


