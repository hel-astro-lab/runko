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
    rhos = np.zeros( (prm.nxfull, prm.ns) )
    
    #sum over velocity dimension
    # wp = sqrt(1/2) ???
    # qn = dv/(q * wp^2 * dx)
    for kk in range(prm.ns):
        for ii in prm.xfull:
            rhos[ii, kk] = np.sum( ff[:, ii, kk] ) * prm.qn[kk]

    #sum over spatial dimension
    rho = np.sum(rhos, 1)

    # XXX why remove the mean?
    rho0 = np.mean( rho[prm.xmid] )
    rho = rho - rho0

    return rho

def wrap(x, prm):

    #left halo to right edge
    x[prm.xLb] = x[prm.xRe]

    #right halo to left edge
    x[prm.xRb] = x[prm.xLe]

    return x


def poisson(ex, rho, prm):

    for ii in prm.xmid:
        ex[ii] = ex[ii-1] + rho[ii]

    ex = wrap(ex, prm)

    ex0 = np.sum( ex[prm.xmid] ) / prm.nx
    ex -= ex0

    #relocate
    #for ii in prm.xmid:
    #    fex[ii] = ( ex[ii] + ex[ii-1] )/2.0

    #fex = wrap(fex, prm)

    return ex



def position_linear(ff, vx, prm):
    ajxs = np.zeros( (prm.nxfull, prm.ns) )
    flux = np.zeros( (prm.nvfull, prm.nxfull) )

    for kk in range(prm.ns):
        aa = vx[:, kk] / prm.dx * prm.dt
        fa = np.floor(aa)

        #XXX fixme
        for ii in prm.xmid:
            iss = np.ones(prm.nvfull)*ii - fa
            flux[:, ii] = aa[:] * np.diag(ff[:, iss.astype(int), kk])

        ff[:, prm.xmid, kk] -= (flux[:, prm.xmid+1] - flux[:, prm.xmid])
                
    #wrap boundaries
    #ff[:, 0:2, :]               = ff[:, prm.nx:prm.nx+2, :]
    #ff[:, prm.nx+3:prm.nx+6, :] = ff[:, 3:6, :]
    ff[:, prm.xLb, :] = ff[:, prm.xRe, :]
    ff[:, prm.xRb, :] = ff[:, prm.xLe, :]
                    
    return ff


def velocity_linear(ff, ex, prm):

    #interpolate half-integer staggered Ex to full integer grid fex
    fex = np.zeros(prm.nxfull)
    for ii in range(1,prm.nxfull):
        fex[ii] = (ex[ii-1] + ex[ii])/2.0
    fex[0] = (ex[-1] + ex[0])/2.0


    flux = np.zeros( (prm.nvfull, prm.nxfull) )

    jj = np.arange(prm.nvfull)
    for kk in range(prm.ns):
        aa = fex[:] * prm.qm[kk]/prm.dv[kk]*prm.dt
        fa = np.floor(aa).astype(int)

        for ii in range(prm.nx+5):
            js = jj - fa[ii] - 1
            flux[jj, ii] = aa[ii] * ff[js, ii, kk]

        ff[1:prm.nv+5, :, kk] = ff[1:prm.nv+5, :, kk] \
                - (flux[1:prm.nv+5, :] - flux[0:prm.nv+4, :])

    return ff


def current(ff, vx, prm):

    ajxs = np.zeros( (prm.nxfull, prm.ns) )
    for kk in range(prm.ns):
        aa = vx[:, kk] / prm.dx * prm.dt

        for ii in prm.xfull:
            ajxs[ii, kk] = np.sum( ff[:, ii, kk] ) * prm.qn[kk]

    ajx = np.sum(ajxs, 1)
    ajx -= np.mean( ajx[prm.xmid] ) 

    return ajx



def efield(ex, ajx, prm):

    #amperes law E_n+1 = E_n - J
    #ex[3:prm.nx+3] = ex[3:prm.nx+3] - ajx[3:prm.nx+3]
    ex[prm.xmid] = ex[prm.xmid] - ajx[prm.xmid]

    #wrap
    #ex[0:2]               = ex[prm.nx:prm.nx+2]
    #ex[prm.nx+3:prm.nx+4] = ex[3:4]
    ex = wrap(ex, prm)

    #shift
    #fex = np.zeros(prm.nxfull)
    #for ii in prm.xmid:
    #    fex[ii] = (ex[ii] + ex[ii-1])/2.0

    #wrap
    #fex[0:2]               = fex[prm.nx:prm.nx+2]
    #fex[prm.nx+3:prm.nx+4] = fex[3:4]
    #fex = wrap(fex, prm)

    return ex


#def efield_f(fex, ajx, prm):
#    #fex[3:prm.nx+3] = fex[3:prm.nx] - ajx[3:prm.nx+3]
#    fex[prm.xmid] = fex[prm.xmid] - ajx[prm.xmid]
#
#    #wrap
#    #fex[0:2]               = fex[prm.nx:prm.nx+2]
#    #fex[prm.nx+3:prm.nx+4] = fex[3:4]
#    fex = wrap(fex, prm)
#
#    return fex




#initialize
#-------------------------------------------------- 
#load configuration
ff, gx, gv, ex, ajx, xx, vx, kx, kv = initial(prm)


#initial step
rho = charge(ff, prm)
ex = poisson(ex, rho, prm)
ff = position_linear(ff, vx, prm)
ajx = current(ff, vx, prm)
ex = efield(ex, ajx, prm)




#-------------------------------------------------- 
# main loop

jtime = 0
for jtime in range(prm.ntime):
    print "-----------", jtime, "----------"
    visz.visualize(fig, jtime, xx, vx, ff, ex, rho, ajx)

    ff  = velocity_linear(ff, ex, prm)
    ff  = position_linear(ff, vx, prm)

    ajx = current(ff, vx, prm)
    rho = charge(ff, prm)

    ex = efield(ex, ajx, prm)
    











