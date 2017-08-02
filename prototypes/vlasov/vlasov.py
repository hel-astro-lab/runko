import numpy as np

from initial import initial
import twostream as prm
import visualize as visz

    
#set up figure
#from pylab import *
import pylab as mlab

fig = mlab.figure(figsize=(10, 12))
mlab.rc('font', family='serif')
mlab.rc('xtick', labelsize='xx-small')
mlab.rc('ytick', labelsize='xx-small')


#--------------------------------------------------
# Random info 
#--------------------------------------------------
#conservative schemes 
# 0 linear
# 1 2nd order
# 2 4th order
# 6 CPIC4
#non-conservative
# 3 - cubic spline
# 5 - CIP3 
#rho = charge(ff, prm)
#ex, fex = poisson(ex, rho, prm)

#3rd possible case ???
#fex = efield_f(fex, ajx, prm)
#ex, fex = efield(ex, ajx, prm)
#--------------------------------------------------


#charge density
def charge(ff, prm):
    rhos = np.zeros( (prm.nx + 6, prm.ns) )
    
    #sum over velocity dimension
    for kk in range(prm.ns):
        for ii in range(prm.nx + 6):
            rhos[ii, kk] = np.sum( ff[:, ii, kk] ) * prm.qn[kk]

    #sum over spatial dimension
    rho = np.sum(rhos, 1)

    # XXX nx+3 or nx+2?
    rho0 = np.mean( rho[3:prm.nx+3] )
    rho = rho - rho0

    return rho

def poisson(ex, rho, prm):

    #XXX +3 or +2
    for ii in range(3, prm.nx+3):
        ex[ii] = ex[ii-1] + rho[ii]

    #wrap boundaries
    ex[0:2]               = ex[prm.nx:prm.nx+2]
    ex[prm.nx+3:prm.nx+4] = ex[3:4]

    ex0 = np.sum( ex[3:prm.nx+3] ) / prm.nx
    ex -= ex0

    #relocate
    fex = np.zeros(prm.nx+6) #already defined
    for ii in range(3, prm.nx+2):
        fex[ii] = ( ex[ii] + ex[ii-1] )/2.0

    fex[0:2]               = fex[prm.nx:prm.nx+2]
    fex[prm.nx+3:prm.nx+6] = fex[3:6]

    return ex, fex

def position_linear(ff, vx, prm):
    ajxs = np.zeros( (prm.nx+6, prm.ns) )
    flux = np.zeros( (prm.nv+6, prm.nx+6) )

    for kk in range(prm.ns):
        aa = vx[:, kk] / prm.dx * prm.dt
        fa = np.floor(aa)

        for ii in range(2, prm.nx+3):
            iss = np.ones(prm.nv+6)*ii - fa
            flux[:, ii] = aa[:] * np.diag(ff[:, iss.astype(int), kk])

        ff[:, 3:prm.nx+3, kk] = ff[:, 3:prm.nx+3, kk] \
                - (flux[:, 3:prm.nx+3] - flux[:, 2:prm.nx+2])
                
    #wrap
    ff[:, 0:2, :]               = ff[:, prm.nx:prm.nx+2, :]
    ff[:, prm.nx+3:prm.nx+6, :] = ff[:, 3:6, :]
                    
    return ff

def velocity_linear(ff, ex, prm):
    flux = np.zeros( (prm.nv+6, prm.nx+6) )

    jj = np.arange(prm.nv+6)
    for kk in range(prm.ns):
        aa = ex[:] * prm.qm[kk]/prm.dv[kk]*prm.dt
        fa = np.floor(aa).astype(int)

        for ii in range(prm.nx+5):
            js = jj - fa[ii] - 1
            flux[jj, ii] = aa[ii] * ff[js, ii, kk]

        ff[1:prm.nv+5, :, kk] = ff[1:prm.nv+5, :, kk] \
                - (flux[1:prm.nv+5, :] - flux[0:prm.nv+4, :])

    return ff

def current(ff, vx, prm):

    ajxs = np.zeros( (prm.nx+6, prm.ns) )
    for kk in range(prm.ns):
        aa = vx[:, kk] / prm.dx * prm.dt

        for ii in range(prm.nx+6):
            ajxs[ii, kk] = np.sum( ff[:, ii, kk] ) * prm.qn[kk]

    ajx = np.sum(ajxs, 1)
    ajx -= np.mean( ajx[3:prm.nx+3] ) # XXX check bounds

    return ajx

def efield(ex, ajx, prm):

    #amperes law E_n+1 = E_n - J
    ex[3:prm.nx+3] = ex[3:prm.nx+3] - ajx[3:prm.nx+3]

    #wrap
    ex[0:2]               = ex[prm.nx:prm.nx+2]
    ex[prm.nx+3:prm.nx+4] = ex[3:4]

    #shift
    fex = np.zeros(prm.nx + 6)
    for ii in range(3, prm.nx+3):
        fex[ii] = (ex[ii] + ex[ii-1])/2.0

    #wrap
    fex[0:2]               = fex[prm.nx:prm.nx+2]
    fex[prm.nx+3:prm.nx+4] = fex[3:4]

    return ex, fex



def efield_f(fex, ajx, prm):
    fex[3:prm.nx+3] = fex[3:prm.nx] - ajx[3:prm.nx+3]

    #wrap
    fex[0:2]               = fex[prm.nx:prm.nx+2]
    fex[prm.nx+3:prm.nx+4] = fex[3:4]

    return fex








#initialize
#-------------------------------------------------- 
#load configuration
ff, gx, gv, ex, fex, ajx, xx, vx, kx, kv = initial(prm)

#initial step
rho = charge(ff, prm)
ex,fex = poisson(ex, rho, prm)
ff = position_linear(ff, vx, prm)
ajx = current(ff, vx, prm)
ex, fex = efield(ex, ajx, prm)




#-------------------------------------------------- 
# main loop

jtime = 0
for jtime in range(prm.ntime):
    print "-----------", jtime, "----------"
    visz.visualize(fig, jtime, xx, vx, ff, fex, rho, ajx)

    ff  = velocity_linear(ff, fex, prm)
    ff  = position_linear(ff, vx, prm)

    ajx = current(ff, vx, prm)
    rho = charge(ff, prm)

    ex, fex = efield(ex, ajx, prm)
    











