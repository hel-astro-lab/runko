import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats
from scipy.optimize import curve_fit

from combine_files import combine_files, combine_tiles
from configSetup import Configuration


#--------------------------------------------------
# combine multiple files
fdir = 'weibel/out/'
conf = Configuration('config-weibel.ini') 

print("files...")
ex   = combine_files(fdir, "field",    "ex",  conf)
#ey   = combine_files(fdir, "field",    "ey",  conf)
#bz   = combine_files(fdir, "field",    "bz",  conf)
rho  = combine_files(fdir, "field",    "rho", conf)
ekin = combine_files(fdir, "analysis", "ekin", conf, isp=0)


#--------------------------------------------------
# read simulation data from file
print "Ex shape:", np.shape(ex)

#Read simulation values
dt = conf.dt*conf.interval
dx = conf.dx

print dt, dx
nt, nx, ny, nz = np.shape(ex)

maxi = -1
time = np.arange(nt)*dt
maxtime = time[maxi]


#--------------------------------------------------
#set up figure
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)

fig = plt.figure(figsize=(3.54, 6.0)) #single column fig
#fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(5, 1, wspace=0.0)


axs = []
axs.append( plt.subplot(gs[0,0]) )
axs.append( plt.subplot(gs[1,0]) )
axs.append( plt.subplot(gs[2,0]) )
axs.append( plt.subplot(gs[3,0]) )
axs.append( plt.subplot(gs[4,0]) )


for ax in axs:
    ax.minorticks_on()
    ax.set_xlabel(r'time $t\omega_{\mathrm{p}}$ ')

    ax.set_xlim((0.0, maxtime))


axs[0].set_ylabel(r'$\ln \delta E_x$')
axs[1].set_ylabel(r'Energy $\epsilon$')
axs[2].set_ylabel(r'$\Delta m$')
axs[3].set_ylabel(r'$\epsilon_K$')
axs[4].set_ylabel(r'$E_T$')


def flatten_spatial(arr):
    return arr.reshape( arr.shape[:-3] + (-1,) )

#--------------------------------------------------
# max{| X |}
axs[0].set_yscale('log')

ex_max = np.max( np.abs( flatten_spatial(ex) ),1 )
axs[0].plot(time, ex_max, 'b-')

#ey_max = np.max( np.abs( flatten_spatial(ey) ),1 )
#axs[0].plot(time, ey_max, 'g-')
#
#bz_max = np.max( np.abs( flatten_spatial(bz) ),1 )
#axs[0].plot(time, bz_max, 'r-')



#axs[0].set_ylim(-20.0, 1.0)


#TODO:
Gm = 0.040
Gms = -22.7 + time*Gm
#axs[0].plot(time, Gms, 'r--')

##################################################
axs[1].set_yscale('log')

ex_edens = np.sum( flatten_spatial(ex*ex), 1 )
axs[1].plot(time, ex_edens)

plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('weibel/weibel.pdf')


#ey_edens = np.sum( flatten_spatial(ey*ey), 1 )
#axs[1].plot(time, ey_edens)
#
#bz_edens = np.sum( flatten_spatial(bz*bz), 1 )
#axs[1].plot(time, bz_edens)


edens = np.sum( flatten_spatial( ex*ex + ey*ey + ez*ez ), 1 )
#edens = np.sum( flatten_spatial( ex*ex ), 1 )
bdens = np.sum( flatten_spatial( bx*bx + by*by + bz*bz ), 1 )
#bdens = np.sum( flatten_spatial( bz*bz ), 1 )

axs[1].plot(time, edens, "b--")
axs[1].plot(time, bdens, "r--")


#TODO: model prediction
#Gms = -16.0 + time*Gm # 1/2 comes from compensation of E_x^2
#axs[1].plot(time, Gms, 'r--')


#axs[1].set_ylim(-10.0, 4.0)

##################################################

prtcls = np.sum( flatten_spatial(rho), 1) #integrate particle density
#prtcls /= prtcls[0]

prtcls = np.abs(prtcls - prtcls[0] )/prtcls[0]
#prtcls = np.clip(prtcls, 1.0e-8, 1.0e2)

axs[2].plot(time, np.log10(prtcls))
#axs[2].plot(time, prtcls)

##################################################

ekintot = np.sum( flatten_spatial(ekin), 1)
axs[3].plot(time, ekintot)


##################################################
print("ekin   max:", np.max(ekintot))
print("efield max:", np.max(edens))
print("bfield max:", np.max(bdens))
print("ratio: ekin/e", np.mean(ekintot)/np.mean(edens))
print("ratio: ekin/b", np.mean(ekintot)/np.mean(bdens))


etot = ekintot + edens + bdens
axs[4].plot(time, etot,    "k-" )
axs[4].plot(time, ekintot, "b--")
axs[4].plot(time, edens,  "r--")
axs[4].plot(time, bdens,  "g--")
#axs[4].plot(time, ex_edens,  "r--")
#axs[4].plot(time, ey_edens,  "r--")
#axs[4].plot(time, bz_edens,  "r--")




plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('weibel/weibel.pdf')
