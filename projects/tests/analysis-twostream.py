import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats



#--------------------------------------------------
# read simulation data from file
f = h5py.File('twostream/out/run.hdf5','r')
#f = h5py.File('twostream/out_ub2/run.hdf5','r')
#f = h5py.File('twostream_kebne/ub1.h5','r')
#f = h5py.File('twostream_kebne/ub4.h5','r')
#f = h5py.File('twostream_kebne/ub4_noclip.h5','r')
#f = h5py.File('twostream_kebne/ub4_hiresv.h5','r')

#f = h5py.File('twostream/out_limits_lowres128/run.hdf5','r')
#f = h5py.File('bump-on-tail/out/run.hdf5','r')

#maxi = 4000

#Read field
ex   = f['fields/Ex'][()]
rho  = f['fields/rho'][()]
ekin = f['fields/ekin'][()]

print("Ex shape:", np.shape(ex))

#Read simulation values
dt = f['params'].attrs['dt']
dx = f['params'].attrs['dx']
print(dt, dx)

nx, ny = np.shape(ex)

time = np.arange(ny)*dt
#maxtime = time[maxi]

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

    #ax.set_xlim((0.0, maxtime))
    #ax.set_xlim((0.0, 50.0))
    #ax.set_xlim((0.0, 70.0))
    #ax.set_xlim((0.0, 54.0))
    #ax.set_xlim((0.0, 400.0))
    #ax.set_xlim((250.0, 400.0))
    #ax.set_xlim((0.0, 1200.0))


axs[0].set_ylabel(r'$\ln \delta E_x$')
axs[1].set_ylabel(r'Energy $\epsilon$')
axs[2].set_ylabel(r'$\Delta m$')
axs[3].set_ylabel(r'$\epsilon_K$')
axs[4].set_ylabel(r'$E_T$')


#time *= 0.321003 #normalize with the growth rate
ex_max = np.max( np.abs(ex),0 )
axs[0].plot(time, np.log(ex_max))

#axs[0].set_ylim(-14.0, 1.0)
#axs[0].set_ylim(-20.0, 1.0)
axs[0].set_ylim(-20.0, -5.0)


#relativistic (vb = 1.0)
#Gm = 0.21
#Gms = -8.0 + time*Gm

#classical (vb = 0.025)
#Gm = 0.24377
#Gms = -8.0 + time*Gm

#Gm = 0.35
#Gms = -18.5 + time*Gm

#Gm = 0.25 (vth 0.1)
#Gm = 0.17
#Gms = -21.0 + time*Gm


#Gm = 0.32 #(vth=0.05)
#Gms = -18.3 + time*Gm 


#bump-on-tail
#Gm = 0.12
#Gm = 0.040
#Gms = -22.7 + time*Gm


#gammab = 1.414
gammab = 2.236
#gammab = 3.162 
#gammab = 4.123

Gi = 1.0/np.sqrt(8.0*gammab**3.0)
#Gms = -14.0 + time*Gi #ub1
Gms = -13.5 + time*Gi #ub2

axs[0].plot(time, Gms, 'r--')

##################################################


wedens = np.log10( np.sum( ex*ex, 0 ) )
axs[1].plot(time, wedens)

#model prediction
#Gm = 0.243771 #vb = 0.025
#Gm = 0.0 + 0.21j

#Gm = 0.21
#Gms = -3.0 + time*Gm*2.0 # 1/2 comes from compensation of E_x^2

Gms = -9.0 + time*Gi # 1/2 comes from compensation of E_x^2
axs[1].plot(time, Gms, 'r--')


#axs[1].set_ylim(-10.0, 4.0)
#axs[1].set_ylim(-18.0, 4.0)

##################################################

prtcls = np.sum(rho, 0) #integrate particle density
#prtcls /= prtcls[0]

prtcls = np.abs(prtcls - prtcls[1] )/prtcls[1]
#prtcls = np.clip(prtcls, 1.0e-8, 1.0e2)

print("prtcls:", prtcls)
axs[2].plot(time, np.log10(prtcls))
axs[2].set_ylim((-8.0, 0.0))
#axs[2].plot(time, prtcls)

##################################################

ekintot = np.sum(ekin, 0)
axs[3].plot(time, np.log10(ekintot))


##################################################
print("ekin   max:", np.max(ekintot))
print("efield max:", np.max(wedens))
print("ratio:", np.mean(ekintot)/np.mean(wedens))

ekintot = ekintot * dx

etot = ekintot + wedens
axs[4].plot(time, np.log10( etot),    "k-" )
axs[4].plot(time, np.log10( ekintot), "b--")
axs[4].plot(time, np.log10( wedens),  "r--")



plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('simulation.pdf')
