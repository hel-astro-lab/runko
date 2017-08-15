import numpy as np

# default values (sometimes modified below)
dx = 0.5 #grid spacing
dt = 0.05 #time step
nx = 100 #spatial grid size
nv = 100 #velocity grid size
ntime = 1000 #time steps
nph = 200 #number of radiation energy bins
nangs = 2 #number of radiation angles




###################################################
# two-stream instability
#ns = 2 #number of species
#dx = 1.0
#dt = 0.1
#wp   = np.array([ 0.707107,  0.707107]) #plasma frequency
#qm   = np.array([-1.0, -1.0]) #charge-to-mass ratio
#vt   = np.array([ 1.0,  1.0]) #thermal velocity
#vd   = np.array([ 2.0, -2.0]) #drift velocity
#vlim = 10.0
#vmax = np.array([ vlim,  vlim])
#vmin = np.array([-vlim ,-vlim])

##################################################
# Buneman instability
#ns = 2 #number of species
#wp   = np.array([ 1.0,  0.1]) #plasma frequency
#qm   = np.array([-1.0,  0.01]) #charge-to-mass ratio
#vt   = np.array([ 1.0,  0.1]) #thermal velocity
#vd   = np.array([ 3.0,  0.0]) #drift velocity
#vmax = np.array([ 10.0, 2.0])
#vmin = np.array([-10,  -2.0])


##################################################
# Harmonic Langmuir waves
ns = 2 #number of species
dx = 1.0
dt = 0.05
nx = 200
nv = 200
ntime = 4000
wp   = np.array([ 1.0,  1.0]) #plasma frequency
qm   = np.array([-1.0, -1.0]) #charge-to-mass ratio
vt   = np.array([ 1.0,  1.0]) #thermal velocity
vd   = np.array([ 0.0,  5.0]) #drift velocity
vmax = np.array([ 10.0, 10.0])
vmin = np.array([-8.0,-8.0])


##################################################
# Electrostatic dispersion relation
#ns = 2 #number of species
#wp   = np.array([ 1.0,  0.2]) #plasma frequency
#qm   = np.array([-1.0,  0.04]) #charge-to-mass ratio
#vt   = np.array([ 1.0,  0.02]) #thermal velocity
#vd   = np.array([ 0.0,  0.0]) #drift velocity
#vmax = np.array([ 10.0, 0.2])
#vmin = np.array([-10.0,-0.2])


##################################################
# Wave propagation
#dx = 0.5
#dt = 0.05
#nx = 400
#nv = 100
#ntime = 4000
#ns = 1 #number of species
#wp   = np.array([ 1.0, ]) #plasma frequency
#qm   = np.array([ 1.0, ]) #charge-to-mass ratio
#vt   = np.array([ 1.0, ]) #thermal velocity
#vd   = np.array([ 0.0, ]) #drift velocity
#vmax = np.array([ 10.0,])
#vmin = np.array([-10.0,])




##################################################
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
# Noise
nmode = 8.0 #number of noise modes
pamp = 0.01 #positive noise amplitude
namp = 0.00 #negative noise amplitude
