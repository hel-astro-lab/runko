import numpy as np
import math
from pylab import *
import scipy
import os, sys

from scipy.special import kv



##################################################
#set up figure
fig = figure(figsize=(12, 6), dpi=200)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(1, 2)
gs.update(hspace = 0.2)
gs.update(wspace = 0.2)

ax1 = subplot(gs[0,0])
ax2 = subplot(gs[0,1])

ax1.set_xlim(3e-2, 1e3)
ax1.set_ylim(1e-4, 3e0)
ax1.set_yscale('log')
ax1.set_xscale('log')

#ax2.set_ylim(1e6, 5e35)
ax2.set_yscale('log')
ax2.set_xscale('log')


# Computational constants
me = 9.1093826e-28  # g
qe = 4.803250e-10   # cgs charge units
c = 2.99792458e10   # cm/s
#r_0 = 2.817940325e-13 # cm
h_pc = 6.6260693e-27        # Planck constant
sigma_T=6.65245873e-25        # Thomson cross-section


R=1e7       # Size of the medium
t_esc=R/c   # Escape (light-crossing) time
tau=1.0     # initial Thomson optical depth
Bfield=1e5  # Magnetic field, [G]




#electron grid
vx = np.linspace(0.1, 15.0, 100)
gamma = np.sqrt( vx**2.0 + 1.0 )
ff = gamma**(-3.0) * vx**2.0



#photon grid
px = np.logspace(-11, 3, 200)
fp = np.zeros(200)

#angle-averaged relativistic synchrotron spectrum
#from Ghisellini et al 1988 
def CS(x, gamma):
    
    #constants and scaling factors
    Bc = me**2*c**3/qe/(h_pc/2.0/pi)
    b = Bfield/Bc
    Ub = Bfield**2/8.0/pi
    #const = 3.0*np.sqrt(3)/pi * (sigma_T/me/c) * Ub * b * h_pc
    const = 3.0*np.sqrt(3)/pi * (sigma_T*c) * Ub * b * h_pc
    xb = x/(3.0*b*gamma**2)


    Kq = kv(4.0/3.0, xb) #one quarter Bessel
    Kt = kv(1.0/3.0, xb) #one thirds Bessel

    return const*xb**2*( Kq*Kt - (3.0/4.0)*xb*(Kq**2 - Kt**2) )




def ph_evolve(fe, vx, fp, px, dt):

    nvx = len(vx)
    nph = len(px)

    xb = qe*Bfield/(2.0*pi*me*c)*h_pc/(me*c**2) # Cyclotron energy

    #logarithmic middle points of the grid 
    lvx = np.log(vx)
    lvx[:-1] += np.diff(lvx)*0.5

    M_ph = np.zeros( (nph, nph) )
    B_ph = np.zeros( (nph) )

    d = 0.5 #make it central scheme

    #determine absorption and emission terms
    for i in range(nph):
        x = px[i]

        jsum = ksum = 0.0
        for j in range(nvx-1):
            z = np.exp( lvx[j] )
            gamma = np.sqrt( z*z + 1.0 )

            dv   = vx[j+1] - vx[j]
            ldv  = np.log(dv)

            fint = (1.0-d)*fe[j+1] + d*fe[j]
            df   = ( fe[j+1] - fe[j] )/ldv

            jsum += dv*fint*CS(x, gamma)
            ksum += dv*CS(x, gamma) * (gamma/z/z) * (3.0*fint - df)


        #emission
        B_ph[i] = jsum/h_pc

        #absorption
        M_ph[i,i] = ksum/x**2


    #self absorption term c*k
    M_ph *= h_pc**2.0/(8.0*pi*(me*c)**3.0)

    #Crank-Nicolson coefficient; now 1/2, i.e., no relaxation
    c_CN = 0.5

    #define matrix of the linear system
    Mf_ph = np.zeros(nph)
    for i in range(nph):
        Mf_ph[:] += M_ph[:, i]*fp[i]
    B_ph += fp/dt - c_CN*Mf_ph 
    
    for k in range(nph):
        M_ph[k,k] = M_ph[k,k]*(1.0 - c_CN) + 1.0/dt

    #then solve
    fpn = linalg.tensorsolve(M_ph, B_ph)


    return fpn





##################################################
dt = 1.0e-4

ax1.plot(vx, ff, linestyle='solid', marker='.')

fp = ph_evolve(ff, vx, fp, px, dt)

#energy flux
Lp = fp * px * 4.0*pi*R**3/(3.0*t_esc)/1.22e6

ax2.plot(px, Lp, linestyle='solid', marker='.')



savefig("simpl_rad.pdf")
