import numpy as np
from pylab import *
from scipy.special import kv




#set up figure
fig = figure(figsize=(6, 4), dpi=200)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(1, 1)
gs.update(hspace = 0.2)
gs.update(wspace = 0.2)

ax1 = subplot(gs[0,0])


me = 9.1093826e-28  # g
qe = 4.803250e-10   # cgs charge units
c = 2.99792458e10   # cm/s
#r_0 = 2.817940325e-13 # cm
h_pc = 6.6260693e-27        # Planck constant
sigma_T=6.65245873e-25        # Thomson cross-section


R=1.0e7       # Size of the medium
t_esc=R/c   # Escape (light-crossing) time
tau=1.0     # initial Thomson optical depth
Bfield=1e5  # Magnetic field, [G]




#angle-averaged relativistic synchrotron spectrum
#from Ghisellini et al 1988 
def CS(x, gamma):
    
    #constants and scaling factors
    Bc = me**2*c**3/qe/(h_pc/2.0/pi)
    b = Bfield/Bc

    Ub = Bfield**2/8.0/pi
    const = 3.0*np.sqrt(3)/pi * (sigma_T/me/c) * Ub * b * h_pc
    #const = 3.0*np.sqrt(3)/pi * (sigma_T*c) * Ub * b * h_pc
    xb = x/(3.0*b*gamma**2)


    Kq = kv(4.0/3.0, xb) #one quarter Bessel
    Kt = kv(1.0/3.0, xb) #one thirds Bessel

    return const*xb**2*( Kq*Kt - (3.0/4.0)*xb*(Kq**2 - Kt**2) )


#angle-averaged relativistic synchrotron spectrum
#from Ghisellini et al 1988 
def CS2(x, gamma):
    
    #constants and scaling factors
    Bc = me**2*c**3/qe/(h_pc/2.0/pi)
    b = Bfield/Bc

    Ub = Bfield**2/8.0/pi
    const = 3.0*np.sqrt(3)/pi * (sigma_T/me/c) * Ub * b * h_pc

    x = x / (3.0*b*gamma**2)

    ya1 = x**(1.0/3.0)
    ya2 = 0.5 * x**(-0.15) * np.exp(-x*2.0)
    ya = 1.0/( (1.0/ya1) + (1.0/ya2) )

    return const*ya




##################################################

ax1.set_xlim(1.0e-11, 1e3)
#ax1.set_ylim(1e-4, 3e0)
ax1.set_ylim(1e-37, 1.0e-34)
ax1.set_yscale('log')
ax1.set_xscale('log')

x = np.logspace(-11, 3, 200)

for gamma in [1.1, 2.0, 5.0, 10.0, 100.0]:
    y = CS(x, gamma)
    ax1.plot(x, y, "b-")

    y2 = CS2(x, gamma)
    ax1.plot(x, y2, "r--")




if False:
    y = CS2(x, 1.0)
    ax1.plot(x, y, "r--")
    
    ya1 = x**(1.0/3.0)
    ax1.plot(x, ya1, "b-")
    
    ya2 = 0.5 * x**(-0.15) * np.exp(-x*2.0)
    #ya2 = 1.0 * x**(0.5) * np.exp(-x)
    ax1.plot(x, ya2, "b-")
    
    ya = 1.0/( (1.0/ya1) + (1.0/ya2) )
    ax1.plot(x, ya, "g-")


savefig("Ghisellini_sync_spec.pdf")
