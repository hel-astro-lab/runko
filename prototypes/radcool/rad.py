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
ax2.set_xlim(1e-11, 1e4)
ax2.set_ylim(1e7, 1e36)


# Computational constants
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



def symlog_grid(x1, x2, xthr, N1, N2):
    vx = np.zeros(N1 + N2)
    vx1 = np.linspace(x1,   xthr, N1)
    vx2 = np.logspace(np.log10(xthr), np.log10(x2),   N2)
    vx[:N1] = vx1
    vx[N1:] = vx2
    return vx


#electron grid
vx = np.linspace(0.1, 100.0, 300)
#vx = symlog_grid(0.1, 1.0e2, 5.0, 50, 50)
#vx = np.logspace(-4, 4, 100)
#vx = np.logspace(-4, 4, 100)

gamma = np.sqrt( vx**2.0 + 1.0 )
ff = gamma**(-3.0) * vx**2.0

#ff = np.zeros( len(vx) )
#ff[70] = 1.0

#normalize
#fz_int=np.trapz(ff, x=vx)
fz_int = 1.0


print "integral:", fz_int
ff *= tau/(fz_int * sigma_T*R)



#photon grid
#x=exp(lnx)
#fx=8.0*pi*(me*c*x)**3.0/(h_pc**3.0*(exp(x/2.0e-5) - 1.0))

px = np.logspace(-11, 3, 200)
fp = np.zeros(200)



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



#Cyclotron energy
# just constants of nature here
def cyclotron_energy():
    return qe*Bfield/(2.0*pi*me*c)*h_pc/(me*c**2) # Cyclotron energy


#transform linear grid to logarithmic; and shift by half
def log_midpoints(arr):
    larr = np.log( arr )
    larr[:-1] += np.diff(larr)*0.5
    return larr


# Compute synchrotron spectra 
def ph_evolve(fe, vx, fp, px, dt):

    nvx = len(vx)
    nph = len(px)

    xb = cyclotron_energy()

    #logarithmic middle points of the grid 
    lvx = log_midpoints(vx)

    M_ph = np.zeros( (nph, nph) )
    B_ph = np.zeros( (nph) )

    d = 1.0

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


#Electron synchrotron cooling
def el_synch_cool(fe, vx, fp, px, dt):

    nvx = len(vx)
    nph = len(px)

    xb = cyclotron_energy()     #cyclotron energy
    bdens = Bfield**2/8.0/np.pi #magnetic energy density

    #Auxiliary stuff for the matrices
    dg = np.zeros(nvx - 1)
    for i in range(nvx-1):
        z     = vx[i]
        zu    = vx[i+1] #originally exp(lnz[i] + d_lnz) = (?) lnz[i+1]
        g     = (z**2  + 1.0)**0.5
        gu    = (zu**2 + 1.0)**0.5
        dg[i] = (zu - z)*(zu + z)/(gu + g)


    #vector of electron momenta in between logarithmic grid points
    zvec = np.exp( log_midpoints(vx) )[:-1]


    #Focker-Planck constants for cooling (equal to \dot\gamma_s)
    lvx  = np.log(vx)
    dlvx = np.diff(lvx) #logarithmic grid cell width
    A_half = np.zeros(nvx - 1)
    A_half = -4.0*sigma_T*bdens*zvec**2 /(3.0*me*c) /dg


    M_el = np.zeros( (nvx, nvx) )
    for i in range( nvx ):
        for i_pr in range(i-1, i+2):
            #tridiagonal terms
            if i == 0:
                if i_pr == (i + 1): M_el[i,i_pr] += A_half[i]
            elif i == nvx-1:
                if i_pr == i:       M_el[i,i_pr] -= A_half[i-1]
            else:
                if i_pr == i:       M_el[i,i_pr] -= A_half[i-1]
                elif i_pr == (i + 1): M_el[i,i_pr] += A_half[i]


    # calculating matrices entering electron equation
    #Mf_el = np.zeros(nvx)
    Mf_el = np.matmul(M_el,fe)

    B_el = fe/dt - c_CN*Mf_el
    M_el = (1.0 - c_CN)*M_el
    for k in range(nvx):
        M_el[k,k] += 1.0/dt

    fzn = np.linalg.tensorsolve(M_el, B_el)
    
    return fzn


# Compute electron cooling coupled to radiative spectra
def el_evolve(fe, vx, fp, px, dt):

    nvx = len(vx)
    nph = len(px)


    #Auxiliary stuff for the matrices
    dg = np.zeros(nvx - 1)
    for i in range(nvx-1):
        z     = vx[i]
        #TODO fix non-uniform grid bug here
        zu    = vx[i+1] #originally exp(lnz[i] + d_lnz) = (?) lnz[i+1]
        g     = (z**2  + 1.0)**0.5
        gu    = (zu**2 + 1.0)**0.5
        dg[i] = (zu - z)*(zu + z)/(gu + g)


    xb = cyclotron_energy()     #cyclotron energy
    bdens = Bfield**2/8.0/np.pi #magnetic energy density

    B_s_diff = np.zeros(nvx)

    lvx = log_midpoints(vx)
    lpx = log_midpoints(px)

    for i in range(nvx-1):
        z = exp(lvx[i])
        gamma = np.sqrt(1.0 + z*z)

        
        #integral over synchrotron emissivity times photon distribution
        Hint = 0.0
        for j in range(0, nph):
            x = px
            Hint += fp[j] * CS(px, gamma)*dlnx/x


    # TODO: implement the rest

    return









##################################################
dt = 1.0e-4

#Crank-Nicolson coefficient; now 1/2, i.e., no relaxation
c_CN = 0.5
#c_CN = 5.0e-6


el_const = sigma_T * R
print "electron normalization: ", el_const
print ff.max()



#plot starting electron distribution
ax1.plot(vx, ff*el_const, linestyle='solid', marker='.')


#step photon distribution
fp = ph_evolve(ff, vx, fp, px, dt)


#step electron distribution
for i in range(20):
    ff = el_synch_cool(ff, vx, fp, px, dt)
#ax1.plot(vx, ff*el_const, color='red', linestyle='dashed', marker='.')
ax1.plot(vx, ff*el_const, color='red', linestyle='dashed')


#energy flux #mc^2 eV to ergs
Lp = fp * px * 4.0*pi*R**3/(3.0*t_esc)/1.22e6
ax2.plot(px, Lp, linestyle='solid', marker='.')


#plot -5/2 powerlaw
xx = np.logspace(-10, -5, 10)
yy = xx**(4.0/3.0)
yy /= yy[0]/1.0e15 

ax2.plot( xx, yy, "k--")




savefig("simpl_rad.pdf")
