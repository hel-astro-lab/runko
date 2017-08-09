import numpy as np

from initial import initial
import conf as prm
from visualize import visualize
import os, sys
from timer import Timer

    
#set up figure
import pylab as mlab

fig = mlab.figure(figsize=(10, 12))
mlab.rc('font', family='serif')
mlab.rc('xtick', labelsize='xx-small')
mlab.rc('ytick', labelsize='xx-small')




#charge density
def charge(ff, prm):
    rhos = np.zeros( (prm.nxfull, prm.ns) )
    
    #sum over velocity dimension
    # qn = dv/(q * wp^2 * dx)
    for kk in range(prm.ns):
        for ii in prm.xfull:
            rhos[ii, kk] = np.sum( ff[:, ii, kk] ) * prm.qn[kk]

    #sum over spatial dimension
    rho = np.sum(rhos, 1)

    return rho


def wrap(x, prm):

    #left halo to right edge
    x[prm.xLb] = x[prm.xRe]

    #right halo to left edge
    x[prm.xRb] = x[prm.xLe]

    return x


def poisson(ex, rho, prm):

    # XXX why remove the mean?
    rho -= np.mean( rho[prm.xmid] )

    for ii in prm.xmid:
        ex[ii+1] = ex[ii] + rho[ii+1]

    ex = wrap(ex, prm)
    ex -= np.sum( ex[prm.xmid] )/ prm.nx

    return ex



def position(ff, vx, ajx, prm):
    ajxs = np.zeros( (prm.nxfull, prm.ns) )
    flux = np.zeros( (prm.nvfull, prm.nxfull) )

    for kk in range(prm.ns):

        aa = vx[:, kk] / prm.dx * prm.dt #compute shift in units of cells
        fa = np.floor(aa) # upwind direction

        for ii in range(2, prm.nx+3):

            #1st order linear upwind flux
            #ss = np.ones(prm.nvfull)*ii - fa
            #iss = ss.astype(int)  # index array needs to be type casted into int before we can use it
            #flux[:, ii] = aa[:] * np.diag(ff[:, iss, kk])

            #second order conservative scheme
            #flux[:, ii] = aa * ( ff[:, ii+1, kk] + ff[:, ii, kk] )*0.5 \
            #      - aa[:]*aa * ( ff[:, ii+1, kk] - ff[:, ii, kk] )*0.5

            #4th order conservative
            flux[:, ii] = aa[:]    * (-ff[:,ii+2,kk]+ 7.0*ff[:,ii+1,kk]+ 7.0*ff[:,ii,kk]-ff[:,ii-1,kk])/12.0 \
                        + aa[:]**2 * ( ff[:,ii+2,kk]-15.0*ff[:,ii+1,kk]+15.0*ff[:,ii,kk]-ff[:,ii-1,kk])/24.0 \
                        + aa[:]**3 * ( ff[:,ii+2,kk]-     ff[:,ii+1,kk]-     ff[:,ii,kk]+ff[:,ii-1,kk])/12.0 \
                        + aa[:]**4 * (-ff[:,ii+2,kk]+ 3.0*ff[:,ii+1,kk]- 3.0*ff[:,ii,kk]+ff[:,ii-1,kk])/24.0


        #add flux as f_i^t+dt = f_i^t - (U_i+1/2 - U_i-1/2)
        ff[:, prm.xmid, kk] -= (flux[:, prm.xmid] - flux[:, prm.xmid-1])

        #flux
        ajxs[prm.xmid, kk] = np.sum( flux[:, prm.xmid], 0) * prm.qn[kk]
                
    #wrap boundaries
    ff[:, prm.xLb, :] = ff[:, prm.xRe, :]
    ff[:, prm.xRb, :] = ff[:, prm.xLe, :]
                    
    #monotonize flux
    np.clip(ff, 0.0, None, out=ff)

    #reduce flux
    ajx[:] = np.sum( ajxs, 1 )
    ajx = wrap(ajx, prm)

    return ff, ajx


def velocity(f, ex, prm):

    #interpolate half-integer staggered Ex to full integer grid fex
    fex = np.zeros(prm.nxfull)
    for ii in prm.xmid:
        fex[ii] = (ex[ii] + ex[ii-1])/2.0
    fex = wrap(fex, prm)

    flux = np.zeros( (prm.nvfull, prm.nxfull) )
    jj = np.arange(2,prm.nv+3)
    for kk in range(prm.ns):
        aa = fex[:] * prm.qm[kk]/prm.dv[kk] * prm.dt #shift in units of phase space cells
        
        #1st order linear upwind scheme
        #fa = np.floor(aa).astype(int) #upwind direction
        #for ii in prm.xfull:
        #for ii in range(2, prm.nx+3):
        #    js = jj - fa[ii]
        #    flux[jj, ii] = aa[ii] * ff[js, ii, kk]

        #2nd order conservative
        #flux[jj, :] = aa[:] * ( ff[jj+1, :, kk] + ff[jj, :, kk] )*0.5 \
        #      - aa[:]*aa[:] * ( ff[jj+1, :, kk] - ff[jj, :, kk] )*0.5

        #4th order conservative 
        flux[jj, :] = aa[:]    * (-ff[jj+2,:,kk]+ 7.0*ff[jj+1,:,kk]+ 7.0*ff[jj,:,kk]-ff[jj-1,:,kk])/12.0 \
                    + aa[:]**2 * ( ff[jj+2,:,kk]-15.0*ff[jj+1,:,kk]+15.0*ff[jj,:,kk]-ff[jj-1,:,kk])/24.0 \
                    + aa[:]**3 * ( ff[jj+2,:,kk]-     ff[jj+1,:,kk]-     ff[jj,:,kk]+ff[jj-1,:,kk])/12.0 \
                    + aa[:]**4 * (-ff[jj+2,:,kk]+ 3.0*ff[jj+1,:,kk]- 3.0*ff[jj,:,kk]+ff[jj-1,:,kk])/24.0

        #add flux as f_i^t+dt = f_i^t - (U_i+1/2 - U_i-1/2)
        ff[1:prm.nv+5, :, kk] -= ( flux[1:prm.nv+5, :] - flux[0:prm.nv+4, :] )

        #monotonize flux
        np.clip(ff, 0.0, None, out=ff)


    return ff



def efield(ex, ajx, prm):

    #remove mean
    ajx -= np.mean( ajx[prm.xmid] )

    #amperes law E_n+1 = E_n - J
    ex[prm.xmid] = ex[prm.xmid] - ajx[prm.xmid]
    ex = wrap(ex, prm)

    return ex




#initialize
#-------------------------------------------------- 
#load configuration
ff, ex, ajx, xx, vx = initial(prm)

#initial step
rho = charge(ff, prm)
ex = poisson(ex, rho, prm)
ff, ajx = position(ff, vx, ajx, prm)
ex = efield(ex, ajx, prm)


#-------------------------------------------------- 
# main loop
visz = visualize("out", xx, vx)
visz.plot(0, ff, ex, ajx, rho) #plot once to create figures

simulation = np.zeros( (prm.nx, prm.ntime+1, 1) )

#-------------------------------------------------- 
#Save to file
import h5py
f = h5py.File("out/run.hdf5", "w")
grp0 = f.create_group("params")
grp0.attrs['dx']    = prm.dx
grp0.attrs['dt']    = prm.dt
grp0.attrs['nx']    = prm.nx
grp0.attrs['nv']    = prm.nv
grp0.attrs['ns']    = prm.ns
grp0.attrs['ntime'] = prm.ntime



grp = f.create_group("fields")
dset_ex = grp.create_dataset("Ex", (prm.nx, prm.ntime+1), dtype='f')



timer = Timer(["total", "lap"])
timer.start("total")

jtime = 0
time = 0.0
for jtime in range(prm.ntime+1):
    #print "-----------", jtime, "/", time, "----------"

    if (jtime % 20 == 0):
        timer.stats("lap")
        visz.plot(jtime, ff, ex, ajx, rho)
        timer.start("lap")

    ff      = velocity(ff, ex, prm)
    ff, ajx = position(ff, vx, ajx, prm)
    rho = charge(ff, prm)
    ex = efield(ex, ajx, prm)
    #ex = poisson(ex, rho, prm)
    
    time += prm.dt

    #simulation[:, jtime, 0] = ex[prm.xmid]
    dset_ex[:,jtime] = ex[prm.xmid]
    
    timer.lap("lap")



timer.stop("total")

timer.stats("total")
