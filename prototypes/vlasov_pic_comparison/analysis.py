import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats
from scipy.optimize import curve_fit


#--------------------------------------------------
# read simulation data from file
f1 = h5py.File('pic11.hdf5','r')
f2 = h5py.File('vlv11.hdf5','r')

#f1 = h5py.File('pic3.hdf5','r')
#f2 = h5py.File('vlv4.hdf5','r')

#f1 = h5py.File('pic10.hdf5','r')
#f2 = h5py.File('vlv10d.hdf5','r')

#Read field
#--------------------------------------------------
class Data:

    ex = None
    rho = None
    ekin = None
    jx = None

    dt = None
    dx = None
    time = None
    x = None


def read(f):

    d = Data()

    #Read field
    d.ex   = f['fields/Ex'][()]
    d.rho  = f['fields/rho'][()]
    d.ekin = f['fields/ekin'][()]
    d.jx   = f['fields/jx'][()]
    
    print "Ex shape:", np.shape(d.ex)
    
    #Read simulation values
    d.dt = f['params'].attrs['dt']
    d.dx = f['params'].attrs['dx']
    print d.dt, d.dx
    
    nx, ny = np.shape(d.ex)
    
    d.x    = np.arange(nx)*d.dx
    d.time = np.arange(ny)*d.dt

    return d

#--------------------------------------------------

pic = read(f1)
vlv = read(f2)


print("PIC----------------")
print(pic.jx)
print("VLV----------------")
print(vlv.jx)
print("RATIO----------------")
print(pic.jx/vlv.jx)
print(np.mean(pic.jx/vlv.jx))

#--------------------------------------------------
#set up figure
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)

fig = plt.figure(figsize=(3.54, 6.0)) #single column fig
#fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(8, 1, wspace=0.0)


axs = []
axs.append( plt.subplot(gs[0,0]) )
axs.append( plt.subplot(gs[1,0]) )
axs.append( plt.subplot(gs[2,0]) )
axs.append( plt.subplot(gs[3,0]) )
axs.append( plt.subplot(gs[4,0]) )
axs.append( plt.subplot(gs[5,0]) )
axs.append( plt.subplot(gs[6,0]) )
axs.append( plt.subplot(gs[7,0]) )


for ax in axs:
    ax.minorticks_on()
    #ax.set_xlabel(r'time $t\omega_{\mathrm{p}}$ ')
    #ax.set_xlim((0.0, maxtime))
    #ax.set_xlim((0.0, 50.0))
    #ax.set_xlim((0.0, 117.0))
    #ax.set_xlim((0.0, 25.0))


#axs[0].set_ylabel(r'$\ln \delta E_x$')
#axs[1].set_ylabel(r'Energy $\epsilon$')
#axs[2].set_ylabel(r'$\Delta m$')
#axs[3].set_ylabel(r'$\epsilon_K$')
#axs[4].set_ylabel(r'$E_T$')


tn = 10
#N = 1.0/0.45
#N = 2.8
N = 1.0

print "N:", N, " -----------"

#axs[0].plot(pic.x, pic.ex[:,tn], "r-")
#axs[0].plot(vlv.x, vlv.ex[:,tn], "b--")
axs[0].plot(vlv.x, vlv.ex[:,tn] / pic.ex[:,tn], "k-")

#axs[1].plot(pic.x, pic.jx[:,tn], "r-")
#axs[1].plot(vlv.x, vlv.jx[:,tn], "b--")

axs[1].plot(vlv.x, vlv.jx[:,tn] / pic.jx[:,tn], "k-")


#--------------------------------------------------
def et(ax, d):
    print "ET-------------"

    A = np.transpose( d.jx )
    (lines, cols) = np.shape(A)
    window = np.hamming(lines).reshape(lines, 1)
    #A *= window
    
    ny, nx = np.shape(A)
    #print "shape:", np.shape(A)
    
    # configuration space parameter
    #x = np.arange(nx)*dx
    x = d.x
    #print "x:"
    #print x
    x1 = 1
    x2 = nx
    
    #temporal guiding vector 
    #t = np.arange(ny)*dt
    t = d.time
    #print "t:"
    #print t
    t1 = 0
    t2 = ny
    
    
    # Change to spectra by considering |F] 
    #F = A*A #energy
    F = A
    print "min/max:", np.min(F), np.max(F)
    
    
    X, T = np.meshgrid(x, t)
    
    #print "nx=",nx
    #print "ny=",ny
    #print np.shape(X)
    #print np.shape(T)
    
    #slow (but more flexible) pcolormesh that takes guiding grid
    #im = ax.pcolormesh(K[w1:w2, k1:k2], W[w1:w2, k1:k2], F[w1:w2, k1:k2],
    #            cmap='plasma', 
    #            #vmin=ff.min(),
    #            #vmax=ff.max(),
    #            )
    
    F = mstats.winsorize(F, limits=[0.01, 0.01])
    vminmax = np.maximum( np.abs(np.min(F)), np.abs(np.max(F)) )
    
    #faster plotting with imshow
    im = ax.imshow(F[t1:t2, x1:x2],
            extent=[x[x1], x[x2-1], t[t1], t[t2-1]],
            origin='lower',
            aspect='auto',
            interpolation='nearest',
            cmap='RdYlGn',
            vmin=-vminmax,
            vmax= vminmax
            )

    return F
    
    
F1 = et(axs[2], pic)    
F2 = et(axs[3], vlv)    

ny, nx = np.shape(F1)
Rat = (F1 / F2)-1.0
Rat = mstats.winsorize(Rat, limits=[0.01, 0.01])
vminmax = np.maximum( np.abs(np.min(Rat)), np.abs(np.max(Rat)) )
vminmax = 0.2
print("vmin/max:", vminmax)

im = axs[4].imshow(Rat[0:ny, 1:nx],
            extent=[vlv.x[1], vlv.x[nx-1], vlv.time[0], vlv.time[ny-1]],
            origin='lower',
            aspect='auto',
            interpolation='nearest',
            cmap='RdYlGn',
            vmin=-vminmax,
            vmax= vminmax
            )
    

#N = 2.4
#N2 = 0.77

#N2 = 1.0
N = 1.0
N2 = 1.0

axs[5].plot(pic.time,    np.sum(    pic.ex*pic.ex, 0), "r-")
axs[5].plot(N2*vlv.time, np.sum(N*N*vlv.ex*vlv.ex, 0), "b--")
#axs[5].set_yscale('log')
    

#N = 3.6
N = 1.0

axs[6].plot(pic.time,    np.sum(    pic.jx*pic.jx, 0), "r-")
axs[6].plot(N2*vlv.time, np.sum(N*N*vlv.jx*vlv.jx, 0), "b--")
    


#N = 0.315
N = 1.0
axs[7].plot(pic.time, np.sum(pic.ekin, 0), "r-")
axs[7].plot(vlv.time, N*np.sum(vlv.ekin, 0), "b--")




plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('comp.pdf')
