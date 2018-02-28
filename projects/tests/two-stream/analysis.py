import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats



#--------------------------------------------------
# read simulation data from file
f = h5py.File('out/run.hdf5','r')

maxi = 1800

#Read field
ex  = f['fields/Ex'][()]
rho = f['fields/rho'][()]
#ex  = f['fields/Ex' ][:,:maxi]
#rho = f['fields/rho'][:,:maxi]

print "Ex shape:", np.shape(ex)

#Read simulation values
dt = f['params'].attrs['dt']
dx = f['params'].attrs['dx']
print dt, dx

nx, ny = np.shape(ex)

time = np.arange(ny)*dt

maxtime = time[maxi]

#ex_max = np.zeros((ny))
#for i in range(ny):
#    exs     = ex[:,i]
#    ex_max[i] = np.max(np.abs(exs))


#--------------------------------------------------
#set up figure
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)

fig = plt.figure(figsize=(3.54, 3.54)) #single column fig
#fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(3, 1, wspace=0.0)


axs = []
axs.append( plt.subplot(gs[0,0]) )
axs.append( plt.subplot(gs[1,0]) )
axs.append( plt.subplot(gs[2,0]) )


for ax in axs:
    ax.minorticks_on()
    ax.set_xlabel(r'time $t\omega_{\mathrm{p}}$ ')

    ax.set_xlim((0.0, maxtime))

axs[0].set_ylabel(r'$\ln \delta E_x$')
axs[1].set_ylabel(r'Energy $\epsilon$')
axs[2].set_ylabel(r'$\Delta m$')


#time *= 0.321003 #normalize with the growth rate
ex_max = np.max( np.abs(ex),0 )
axs[0].plot(time, np.log(ex_max))

##################################################

wedens = np.log10( np.sum( ex*ex, 0 ) )
axs[1].plot(time, wedens)

##################################################

prtcls = np.sum(rho, 0)*dx #integrate particle density
#prtcls /= prtcls[0]


prtcls = np.abs(prtcls - prtcls[0] )/prtcls[0]
prtcls = np.clip(prtcls, 1.0e-8, 1.0e2)

axs[2].plot(time, np.log10(prtcls))
#axs[2].plot(time, prtcls)



plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('twostream.pdf')
