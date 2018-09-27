import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats

from combine_files import combine_files, combine_tiles
from configSetup import Configuration

def flatten_spatial(arr):
    return arr.reshape( arr.shape[:-3] + (-1,) ).T



read_files = True
read_runf = False

#--------------------------------------------------
# read simulation data from file
#fdir = 'twostream/out/'
fdir = 'twostream_kebne/ub2/'


if read_files:
    #conf = Configuration('config-twostream-relativistic.ini') 
    conf = Configuration(fdir+'config-ts2.ini') 
    #--------------------------------------------------
    
    print("files...")
    ex_sparse =  flatten_spatial( combine_files(fdir, "field",    "ex",   conf) ) 
    #ey_sparse =  flatten_spatial( combine_files(fdir, "field",    "ey",   conf) )
    #ez_sparse =  flatten_spatial( combine_files(fdir, "field",    "ez",   conf) )
    #bx_sparse =  flatten_spatial( combine_files(fdir, "field",    "bx",   conf) )
    #by_sparse =  flatten_spatial( combine_files(fdir, "field",    "by",   conf) )
    #bz_sparse =  flatten_spatial( combine_files(fdir, "field",    "bz",   conf) )
    rh_sparse =  flatten_spatial( combine_files(fdir, "field",    "rho",  conf) )
    ed_sparse1 = flatten_spatial( combine_files(fdir, "analysis", "edens",conf, isp=0) )
    ed_sparse2 = flatten_spatial( combine_files(fdir, "analysis", "edens",conf, isp=1) )
    ed_sparse = ed_sparse1 + ed_sparse2

    dt = conf.dt
    dx = conf.dx

    nx, nt = np.shape(ex_sparse)
    time_s = np.arange(nt)*dt*conf.interval
    
    print("nt,nx", np.shape(ex_sparse))
    print("time_s:", time_s)


if read_runf:
    f = h5py.File(fdir+'run.hdf5','r')
    dt = f['params'].attrs['dt']
    dx = f['params'].attrs['dx']
    print(dt, dx)


    #--------------------------------------------------
    #Read field
    ex   = f['fields/Ex'][()]
    rho  = f['fields/rho'][()]
    ekin = f['fields/ekin'][()]
    
    print("Ex shape:", np.shape(ex))
    
    #Read simulation values
    nx, ny = np.shape(ex)
    time = np.arange(ny)*dt




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
    ax.set_xlim((0.0, 50.0))


axs[0].set_ylabel(r'$\ln \delta E_x$')
axs[1].set_ylabel(r'Energy $\epsilon$')
axs[2].set_ylabel(r'$\Delta m$')
axs[3].set_ylabel(r'$\epsilon_K$')
axs[4].set_ylabel(r'$E_T$')


if read_runf:

    #time *= 0.321003 #normalize with the growth rate
    ex_max = np.max( np.abs(ex),0 )
    axs[0].plot(time, np.log(ex_max))
    
    
    

if read_files:
    ex_max_s = np.max( np.abs(ex_sparse), 0)
    axs[0].plot(time_s, np.log(ex_max_s))


##################################################
# growth rates

#gammab = 1.414
gammab = 2.236
#gammab = 3.162 
#gammab = 4.123
    
Gi = 1.0/np.sqrt(8.0*gammab**3.0)
#Gms = -14.0 + time*Gi #ub1

if read_runf:
    Gms = -8.5 + time*Gi #ub2
    axs[0].plot(time, Gms, 'r--')

if read_files:
    Gms = -8.5 + time_s*Gi #ub2
    axs[0].plot(time_s, Gms, 'r--')


#axs[0].set_ylim(-14.0, 1.0)
#axs[0].set_ylim(-20.0, 1.0)
axs[0].set_ylim(-10.0, -5.0)


##################################################

if read_runf:
    wedens = np.log10( np.sum( ex*ex, 0 ) )
    axs[1].plot(time, wedens)
    
    #model prediction
    #Gm = 0.243771 #vb = 0.025
    #Gm = 0.0 + 0.21j
    
    #Gm = 0.21
    #Gms = -3.0 + time*Gm*2.0 # 1/2 comes from compensation of E_x^2
    
    Gms = -5.2 + time*Gi # 1/2 comes from compensation of E_x^2
    axs[1].plot(time, Gms, 'r--')
    
    #axs[1].set_ylim(-10.0, 4.0)
    #axs[1].set_ylim(-18.0, 4.0)


if read_files:
    wedens_s = np.sum( ex_sparse**2.0, 0)
    axs[1].plot(time_s, np.log10(wedens_s))

    Gms = -5.2 + time_s*Gi # 1/2 comes from compensation of E_x^2
    axs[1].plot(time_s, Gms, 'r--')


##################################################
axs[2].set_yscale('log')

if read_runf:
    prtcls = np.sum(rho, 0) #integrate particle density
    #prtcls /= prtcls[0]
    
    prtcls = np.abs(prtcls - prtcls[1] )/prtcls[1]
    #prtcls = np.clip(prtcls, 1.0e-8, 1.0e2)
    
    print("prtcls:", prtcls)
    axs[2].plot(time, prtcls)
    axs[2].set_ylim((-8.0, 0.0))
    #axs[2].plot(time, prtcls)

if read_files:
    dens = np.sum(rh_sparse, 0)
    dens = np.abs(dens - dens[1])/dens[1] #normalize to initial value

    axs[2].plot(time_s, dens)
    axs[2].set_ylim((1.0e-7, 1.0e1))



##################################################
# units
ef_units = dt*2.
ef_units = 8.0*dx*(dt/dx)**2.
print("ef_units", ef_units)


if read_runf:
    ekintot = np.sum(ekin, 0)
    axs[3].plot(time, np.log10(ekintot))


    ##################################################
    print("ekin   max:", np.max(ekintot))
    print("efield max:", np.max(wedens))
    print("ratio:", np.mean(ekintot)/np.mean(wedens))
    
    ekin = np.sum(ekin, 0) #E_kin density
    
    eef = np.sum(ex*ex, 0) #E_x density
    eef /= ef_units
    
    etot = ekin + eef
    axs[4].plot(time, etot,    "k-" )
    axs[4].plot(time, ekin, "b--")
    axs[4].plot(time, eef,  "r--")


if read_files:
    # sparse sampling
    ekin_s = np.sum(ed_sparse,0)
    eef_s  = np.sum(ex_sparse**2.0,0)/ef_units
    etot_s = ekin_s + eef_s
    
    axs[4].plot(time_s, etot_s, "k." )
    axs[4].plot(time_s, ekin_s, "b.")
    axs[4].plot(time_s, eef_s,  "r.")



plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('simulation.pdf')
