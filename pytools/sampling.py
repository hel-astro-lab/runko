# -*- coding: utf-8 -*- 

import numpy as np


# 2D velocity
#Change from u = |u| to (ux, uy, 0)
def velxy(u):

    #now we have valid u = abs(u_i)
    #x1 = np.random.rand()
    x2 = np.random.rand()
    #x3 = np.random.rand()

    #2d treatment; collapsing ux and flipping ux <-> uz
    ux = u*np.cos(2.0*np.pi*x2)
    uy = u*np.sin(2.0*np.pi*x2)
    uz = 0.0

    return ux, uy, uz

# 3D velocity
#Change from u = |u| to (ux, uy, uz)
def velxyz(u):

    #now we have valid u = abs(u_i)
    x1 = np.random.rand()
    x2 = np.random.rand()
    x3 = np.random.rand()

    #3d treatment
    ux = u*(2.0*x1 -1.0)
    uy = 2.0*u*np.sqrt(x1*(1.0-x1))*np.cos(2.0*np.pi*x2)
    uz = 2.0*u*np.sqrt(x1*(1.0-x1))*np.sin(2.0*np.pi*x2)

    return ux, uy, uz



#Sobol method to sample 4-velocities u = g*v
def sobol_method(T):
    x4 = np.random.rand()
    x5 = np.random.rand()
    x6 = np.random.rand()
    x7 = np.random.rand()

    u = -T*np.log(x4*x5*x6)
    n = -T*np.log(x4*x5*x6*x7)

    if n*n - u*u < 1.0:
        return sobol_method(T)

    return u

#Box-Muller sampling
def BoxMuller_method(T):
    vth = np.sqrt(2.0*T)

    #Box-Muller sampling
    rr1 = np.random.rand()
    #rr2 = np.random.rand()
    #r1  = np.sqrt(-2.0*np.log(rr1))*np.cos(2.0*np.pi*rr2)
    #r2  = np.sqrt(-2.0*np.log(rr1))*np.sin(2.0*np.pi*rr2)
    #ux = r1*vth 
    #uy = r2*vth

    return np.sqrt(-2.0*np.log(rr1))*vth


# Boosted Maxwellian according to Zenitani 2015
# theta = kT/mc^2
# Gamma = Lorentz boosting factor (or if G<1 it is interpreted as beta)
# direction = -1/+1 for -/+x; -2/+2 for -/+y; -3/+3 for -/+z
# dimensionality = 2 for xy; 3 for xyz
def sample_boosted_maxwellian(theta, Gamma, direction=1, dims=2):

    #For relativistic case we use Sobol method, inverse method otherwise
    if theta > 0.2:
        u = sobol_method(theta)
    else:
        u = BoxMuller_method(theta)


    #now get components
    if dims == 2:
        ux, uy, uz = velxy(u)
    if dims == 3:
        ux, uy, uz = velxyz(u)


    #if no drift, we just return the non-boosted distribution
    if Gamma == 0:
        return ux, uy, uz, u

    if Gamma < 1.0:
        #We interpret this as v/c = beta
        beta = Gamma
        Gamma = 1.0/np.sqrt(1.0 - beta*beta)
    else: 
        #else as bulk lorentz factor
        beta = np.sqrt(1.0 - 1.0/Gamma/Gamma)

    #beta = 1.0/sqrt(1.0 + Gamma*Gamma)
    X8 = np.random.rand()
    vx = ux/np.sqrt(1.0 + u*u)
    if -beta*vx > X8:
        ux = -ux

    ux = Gamma*(ux + beta*np.sqrt(1.0 + u*u))
    u = np.sqrt(ux*ux + uy*uy + uz*uz)


    # swap directions
    if   direction == -1:
        ux = -ux
    elif direction == +1:
        ux = ux
    elif direction == -2:
        tmp = -ux
        ux = uy
        uy = tmp
    elif direction == +2:
        tmp = +ux
        ux = uy
        uy = tmp
    elif direction == -3:
        tmp = -ux
        ux = uz
        uz = tmp
    elif direction == +3:
        tmp = +ux
        ux = uz
        uz = tmp
    else:
        raise Exception('Invalid direction; give |d| <= 3')

    return ux, uy, uz, u


def sample_blackbody(delgam):

    # TODO
    uu = 1.0
    ux,uy,uz = velxyz(uu)

    return ux, uy, uz, uu



##################################################

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    plt.fig = plt.figure(1, figsize=(4,3))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    for ai in range(1):
        axs.append( plt.subplot(gs[ai]) )


    #testing sobol vs box-muller
    if False:
        Gamma = 0.0
        #T = 2.0e-4
        T = 0.2

        N = 10000
        n1 = np.zeros(N)
        n2 = np.zeros(N)
        n3 = np.zeros(N)
        for n in range(N):
        
            #n1[n] = rejection_sampling(A)
        
            #Sobol for relativistic
            #vx, vy, vz, u = sobol_method(T)
            #n2[n] = u
        
        
            #non rel maxwell with rejection method
            #vx, vy, vz = thermal_plasma(T)
            #u = sqrt(vx*vx + vy*vy + vz*vz)
            #n3[n] = u
        
            #n1[n] = drifting_rejection(A, drift)
        
            ux, uy, uz, u = sample_boosted_maxwellian(T, Gamma, direction=-1)
            #p = u/sqrt(1.0-u*u)
            n2[n] = uy
        
            ux, uy, uz, u = sample_boosted_maxwellian(T+0.0001, Gamma, direction=-1)
            #p = u/sqrt(1.0-u*u)
            n1[n] = uy
        
        #plot
        #print n2
        
        axs[0].hist(n1, 100, color="black", alpha=0.3, density=True)
        axs[0].hist(n2, 100, color="red"  , alpha=0.3, density=True)
    

    #testing Sobol with relativistic temperatures
    if True:
        from scipy.special import kn

        Gamma = 0.0
        #T = 2.0e-4
        T = 0.3

        N = 10000
        n1 = np.zeros(N)
        n2 = np.zeros(N)
        n3 = np.zeros(N)
        for n in range(N):
            ux, uy, uz, u = sample_boosted_maxwellian(T, Gamma, direction=-1)
            gamma = np.sqrt(1.0 + ux*ux + uy*uy + uz*uz)
            n1[n] = gamma
        
            ux, uy, uz, u = sample_boosted_maxwellian(T+0.0001, Gamma, direction=-1)
            gamma = np.sqrt(1.0 + ux*ux + uy*uy + uz*uz)
            n2[n] = gamma
        

        axs[0].set_xscale('log')
        #axs[0].set_yscale('log')
        axs[0].set_xlim((1.0, 10.0))

        axs[0].hist(n1, 100, color="black", alpha=0.3, density=True)
        axs[0].hist(n2, 100, color="red", alpha=0.3, density=True)
    
        #gamma grid (with beta)
        gs = np.logspace(0.0, 2.0, 100) + 0.01
        beta = np.sqrt(1.0 - 1.0/gs**2)

        #maxwell-juttner distribution
        K2T = kn(2, 1.0/T) #modified Bessel function of the second kind (with integer order 2)
        mwj = (gs**2 * beta)/(T*K2T)*np.exp(-gs/T)

        axs[0].plot(gs, mwj, 'r-')
    
    
    fname = 'maxwells.pdf'
    plt.savefig(fname)





