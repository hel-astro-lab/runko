import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats



#--------------------------------------------------
# read simulation data from file
f = h5py.File('out/run.hdf5','r')


#Read field
ex_f = f['fields/Ex']
ex_f = ex_f[:,50:]
print "Ex shape:", np.shape(ex_f)

#Read simulation values
dt = f['params'].attrs['dt']
dx = f['params'].attrs['dx']
print dt, dx

nx, ny = np.shape(ex_f)

time = np.arange(ny)*dt

lnex = np.zeros((ny))
for i in range(ny):
    exs     = ex_f[:,i]
    lnex[i] = np.log10( np.max(exs) )


#--------------------------------------------------
#set up figure
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)

fig = plt.figure(figsize=(3.54, 3.54)) #single column fig
#fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(2, 1, wspace=0.0)


axs = []
axs.append( plt.subplot(gs[0,0]) )
axs.append( plt.subplot(gs[1,0]) )


for ax in axs:
    ax.minorticks_on()
    ax.set_xlabel(r'time $t\omega_{\mathrm{p}}$ ')

axs[0].set_ylabel(r'$\log_{10} \delta E_x$')
axs[1].set_ylabel(r'Energy $\epsilon$')

axs[0].plot(time, lnex)


##################################################

wedens = np.log10( np.sum( ex_f*ex_f, 0 ) )
axs[1].plot(time, wedens)



plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('landau.pdf')
