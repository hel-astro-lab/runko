import numpy as np
from pylab import *
from scipy.special import kv


from scipy import interpolate



#set up figure
fig = figure(figsize=(6, 4), dpi=200)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(1, 2)
gs.update(hspace = 0.2)
gs.update(wspace = 0.2)

ax1 = subplot(gs[0,0])
ax2 = subplot(gs[0,1])


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
def CSapprox(x, gamma):
    
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


def el_cool(vx, ff):

    vxn = np.zeros(len(vx))
    ffn = np.zeros(len(ff))
    dv = vx[2] - vx[1]
    dt = 0.01
    for i, v in enumerate(vx):
        gamma = np.sqrt(1.0 + v*v)
        beta  = v / gamma

        dgamma = gamma**2 * beta**2
        gamman = gamma - dgamma*dt
        vn = np.sqrt(gamman**2 - 1.0)

        #fa = np.int( vn/dv )

        vxn[i] = vn

    ffn = ff
    return vxn, ffn


def el_cool2(vx, ff):

    vxn = np.zeros(len(vx))
    ffn = np.zeros(len(ff))
    dv = vx[2] - vx[1]
    dt = 0.02
    for i, v in enumerate(vx):
        gamma = np.sqrt(1.0 + v*v)
        beta  = v / gamma

        dgamma = gamma**2 * beta**2
        gamman = gamma - dgamma*dt
        vn = np.sqrt(gamman**2 - 1.0)

        vxn[i] = vn


    lff  = np.log10(ff)
    f = interpolate.interp1d(vxn, lff, fill_value='extrapolate')
    fnew = f(vx)
    return vx, 10.0**fnew


    #re-grid to the old momenta grid
    lvxn = np.log10(vxn)
    lvx  = np.log10(vx)
    lff  = np.log10(ff)

    f = interpolate.interp1d(lvxn, lff, fill_value='extrapolate')
    fnew = f(lvx)

    return vx, 10.0**fnew



##################################################

ax1.set_xlim(1.0e-1, 3e1)
#ax1.set_ylim(1e-4, 3e0)
#ax1.set_ylim(1e-37, 1.0e-34)
ax1.set_yscale('log')
ax1.set_xscale('log')

#ax1.set_xlim(1.0e-1, 3e1)
#ax1.set_ylim(1e-4, 3e0)
#ax1.set_ylim(1e-37, 1.0e-34)
ax2.set_yscale('log')
ax2.set_xscale('log')

##################################################



vx = np.linspace(0.1, 20.0, 50)
gamma = np.sqrt( vx**2.0 + 1.0 )
ff = gamma**(-3.0) * vx**2.0

ax1.plot(vx, ff, "k-", marker=".")
ax2.plot(gamma, ff, "k-", marker=".")



vx2, ff2 = el_cool(vx, ff)

gamma2 = np.sqrt(vx2**2 + 1.0)
ax1.plot(vx2, ff2, "r--")
ax2.plot(gamma2, ff2, "r-")


# re-gridding
vx3, ff3 = el_cool2(vx, ff)
gamma3 = np.sqrt(vx3**2 + 1.0)
ax1.plot(vx3, ff3, "b-", alpha=0.7)
ax2.plot(gamma3, ff3, "b-", alpha=0.7)



savefig("sync_cooling.pdf")
