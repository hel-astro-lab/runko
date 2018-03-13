import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats
from scipy.optimize import curve_fit


#--------------------------------------------------
# read simulation data from file
f = h5py.File('out/run.hdf5','r')
#f = h5py.File('out_khat045/run.hdf5','r')
#f = h5py.File('out_khat014/run.hdf5','r')

maxi = 100

#Read field
ex   = f['fields/Ex'][()]
rho  = f['fields/rho'][()]
ekin = f['fields/ekin'][()]
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
    #ax.set_xlim((0.0, 117.0))
    ax.set_xlim((0.0, 50.0))


axs[0].set_ylabel(r'$\ln \delta E_x$')
axs[1].set_ylabel(r'Energy $\epsilon$')
axs[2].set_ylabel(r'$\Delta m$')
axs[3].set_ylabel(r'$\epsilon_K$')
axs[4].set_ylabel(r'$E_T$')


#time *= 0.321003 #normalize with the growth rate
ex_max = np.max( np.abs(ex),0 )
axs[0].plot(time, np.log10( ex_max ))
#axs[0].set_ylim((-20, 10))





##################################################

wedens = np.sum( ex*ex, 0 ) #*2.0*np.pi
axs[1].plot(time, np.log10(wedens))


#measure frequency
#w = np.fft.fft(ex_max)
#freqs = np.fft.fftfreq(len(time))
#
#idx = np.argmax(np.abs(w))
#freq = freqs[idx]
#print("max frq:", freq)
#
#axs[1].plot(freqs, np.log10(w*w) )



#plot linear analysis results
#omega=1.415661888604536 - 0.153359466909605j
#omega=1.415 - 0.153j

#parameters
delgam = 0.001
Nw = 2.0
khat = Nw * 2.0*np.pi*np.sqrt(delgam)/(nx*dx)
khat = 0.14

#real part
wr = 1.0 + ((3.0/2.0)*khat**2 + (15.0/8.0)*khat**4 + (147/16.0)*khat**6
           + 736.437*khat**8
           - 14729.3*khat**10
           + 105429*khat**12.0
           - 370151*khat**14.0
           + 645538*khat**16.0
           - 448190*khat**18.0)

#imaginary part
wi = -0.5*np.sqrt(np.pi/2.0)*(
        1.0/khat**3 - 6.0*khat 
        - 40.7173*khat**3.0
        - 3900.23*khat**5.0
        - 2462.25*khat**7.0
        - 274.99*khat**9.0)
wi *= np.exp(-0.5/khat**2 - 3.0/2.0 
        - 3.0*khat**2.0
        - 12.0*khat**4.0
        - 575.51*khat**6.0
        + 3790.16*khat**8.0
        - 8827.54*khat**10.0
        + 7266.87*khat**12.0)
#omega= wr + wi*1.0j

print("khat:",khat)
print("wr:",wr)
print("wi:",wi)


##################################################
omega = 1.3 - 0.10j #khat=0.45
#omega = 1.0004 #khat = 0.14
#omega = 0.99

def landau_osc(t, omega):
    f = 0.5*wedens[0]*np.real( np.exp(-1j*omega*(t-0.3)))**2.0
    return f

#popt, pcov = curve_fit(landau_osc, time, wedens, p0=omega, bounds=(0.5, 2.0))
#omega = popt
#print(popt)



#line 
tskip = 0.28 #0.8
lin_analysis  = wedens[0]*np.abs(  np.exp(-1j*omega*(time-tskip)))**2.0

#with frequency
lin_analysis2 = wedens[0]*np.real( np.exp(-1j*omega*(time-tskip)))**2.0


#% linear analysis
#plot(time,0.5*fieldenergy(1)*abs(exp(-1j*omega*(time-0.4))).^2);
#% linear analysis with frequency
#plot(time,0.5*fieldenergy(1)*real(exp(-1j*omega*(time-0.4))).^2);

axs[1].plot(time, np.log10( lin_analysis  ), 'g-', alpha=0.6)
axs[1].plot(time, np.log10( lin_analysis2 ), 'r-', alpha=0.6)


##################################################

prtcls = np.sum(rho, 0) #integrate particle density
#prtcls /= prtcls[0]


prtcls = np.abs(prtcls - prtcls[0] )/prtcls[0]
#prtcls = np.clip(prtcls, 1.0e-8, 1.0e2)

axs[2].plot(time, np.log10(prtcls))
#axs[2].plot(time, prtcls)


##################################################

ekintot = np.sum(ekin, 0) #*n0
axs[3].plot(time, np.log10(ekintot))
#axs[3].plot(time, wedens, 'r-')


##################################################
#etot = np.sum(ex*ex/(8.0*np.pi) ,0) #+ np.sum(ekin,0)
#etot = np.sum(ekin,0)/np.sum(rho, 0)
#etot = np.sum(ekin,0)*dx
#etot = etot/etot[0]

print("ekin   max:", np.max(ekintot))
print("efield max:", np.max(wedens))
print("ratio:", np.mean(ekintot)/np.mean(wedens))


#etot = ekintot + wedens
#etot /= etot[0]
#etot = np.abs(etot/etot[0] - 1.0)
#axs[4].plot(time, np.log10( etot) )

#ekintot = ekintot * dx*dx #try normalization

ekintot = ekintot * dx

etot = ekintot + wedens
axs[4].plot(time, np.log10( etot),    "k-" )
axs[4].plot(time, np.log10( ekintot), "b--")
axs[4].plot(time, np.log10( wedens),  "r--")


plt.subplots_adjust(left=0.18, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('simulation.pdf')
