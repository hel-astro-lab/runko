import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy.stats import special.erf as erf



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
    ax.set_xlabel(r'x')
    #ax.set_xlim((0.0, 100.0))

axs[0].set_ylabel(r'y')


#complex error function
def Erfi(v):
    return -1.0j*erf(v)


def ts_disp(khat, zb, omega):
    """ two-stream instability dispersion relation Shalaby et al. 2017 """
    vp = omega/khat

    ret = np.sqrt(np.pi/8.0)*(
    (vp + zb)*( Erfi( (vp+zb)/np.sqrt(2) ) -1.0j ) * np.exp(-vp*zb) +
    (vp - zb)*( Erfi( (vp-zb)/np.sqrt(2) ) -1.0j ) * np.exp(+vp*zb)
    )*np.exp(-(vp**2.0 + zb**2.0)/2.0) )

    return np.sqrt(ret + 1.0)



def dispersion_relation(vb, delgam):
    """ solves dispersion relation for range of k:s """

    omega = 1.4 - 0.01j      #initial guess
    zb = vb/np.sqrt(delgam)  #thermal velocity

    karr = np.linspace(2.0, 30.0, 10)
    for i,k in enumerate(karr):
        khat = k*np.sqrt(delgam)




xarr = np.linspace(0.1, 25.0, 20)
yarr = np.zeros((len(xarr))
for i, x in enumerate(xarr):
    yarr[i] = ts_disp


axs[0].plot(xarr, yarr, 'r-')




