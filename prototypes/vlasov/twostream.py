import numpy as np

dx = 1.0
dt = 0.1
nx = 200
nv = 100
ntime = 2000

#species
ns = 2
wp   = np.array([ 0.707107,  0.707107  ])
qm   = np.array([-1.000000, -1.000000])
vt   = np.array([ 1.000000,  1.000000, ])
vd   = np.array([ 0.000000,  4.000000, ])
vmax = np.array([ 9.999990,  9.999900, ])
vmin = np.array([-5.000000, -5.000000])


#originally in the initial() function
dv = np.zeros(ns)
qn = np.zeros(ns)


#noise  3 - white noise
#noise = 3 


nmode = 1.000000;
pamp = 0.01
namp = 0.00
#pphs = 0.00
#nphs = 0.00


#diagnostics
#nplot    = 1000.000000
#emax     = 1.000000
#diagtype = [10.000000, 11.000000, 2.000000, 9.000000, 5.000000, 1.000000 ]


