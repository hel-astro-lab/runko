import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats



#--------------------------------------------------
# read simulation data from file
f = h5py.File('out/run.hdf5','r')

maxi = 100

#Read field
ex   = f['fields/Ex'][()]
rho  = f['fields/rho'][()]
ekin = f['fields/ekin'][()]

print "Ex shape:", np.shape(ex)

#Read simulation values
dt = f['params'].attrs['dt']
dx = f['params'].attrs['dx']
print dt, dx

nx, ny = np.shape(ex)

time = np.arange(ny)*dt
xarr = np.arange(nx)*dx

maxtime = time[maxi]


#--------------------------------------------------
#set up figure
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)

fig = plt.figure(figsize=(3.54, 6.0)) #single column fig
#fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(1, 1, wspace=0.0)


axs = []
axs.append( plt.subplot(gs[0,0]) )

for ax in axs:
    ax.minorticks_on()
    ax.set_xlabel(r'time $t\omega_{\mathrm{p}}$ ')
    #ax.set_xlim((0.0, maxtime))


axs[0].set_ylabel(r'$E_x$ Error')


################################################### 
# calculate field error
#axs[0].set_ylim((-20, 10))

dex = np.zeros((nx, ny))

#dex[:-1,:] = -(ex[1:, :]  - ex[:-1, :])/dx #asymm
#dex[1:,:]  = (ex[:-1, :] - ex[1:, :])/dx
#dex[:-2,:] = (ex[:-2, :] - ex[0:-2, :])/dx

dex[:-1,:]  = (ex[:-1, :] - ex[1:, :])/dx #nabla+
#dex[1:,:]  = (ex[1:, :] - ex[:-1, :])/dx #nabla-


#dex[-1, :] = (ex[-1,:] - ex[0,:])/dx #periodic boundary 



#rhomod = (rho - 1.0)/0.45
rhomod = (rho - 1.0)

it = 0
print( np.max(rhomod[:,it]) / np.max(dex[:,it]))

axs[0].plot(xarr, dex[:,it], 'b-')
axs[0].plot(xarr, rhomod[:,it], 'r-')




plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('charge.pdf')
