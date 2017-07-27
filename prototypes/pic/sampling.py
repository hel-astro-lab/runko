import numpy as np
import math
from pylab import *
import os, sys
import scipy
from scipy.stats import gaussian_kde
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import kn

#set seed to get reproducible errors & results
np.random.seed(0)



#set up figure
fig = figure(figsize=(10, 12), dpi=200)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')


gs = GridSpec(1, 1)
gs.update(hspace = 0.2)
gs.update(wspace = 0.2)

ax1 = subplot(gs[0,0])
#ax1.set_xlim((0, 5))
ax1.set_xlabel(r'velocity $v$')



#def thermal_plasma(theta):
#    fmax = 1.0
#    vmin = -5.0*theta
#    vmax = 5.0*theta
#    
#    vf = vmin + (vmax-vmin)*np.random.rand()
#    f = 0.5*(exp(-(vf*vf)/(2.0*theta*theta)))
#
#    x = fmax*np.random.rand()
#
#    if x > f:
#        return thermal_plasma(theta)
#
#    #now we have valid u = abs(u_i)
#    x1 = np.random.rand()
#    x2 = np.random.rand()
#    #x3 = np.random.rand()
#
#    vx = vf*(2*x1 -1)
#    vy = 2*vf*sqrt(x1*(1-x1))
#
#    #3d treatment
#    #vy = 2*u*sqrt(x1*(1-x1))*cos(2*pi*x2)
#    #vz = 2*u*sqrt(x1*(1-x1))*sin(2*pi*x2)
#    vz = 0.0
#
#    return vx,vy,vz
#
#
##relativistic maxwell-Juttner distribution
## theta is dimensionless temperature
#def thermal_rel_plasma(theta):
#
#    fmax = 1.0/kn(2,1.0/theta)
#    vmin = -20.0*theta
#    vmax = 20.0*theta
#    vf = vmin + (vmax-vmin)*np.random.rand()
#    
#    f = exp(-sqrt(1+vf*vf)/theta)*vf*vf
#
#    x = fmax*np.random.rand()
#
#    if x > f:
#        return thermal_rel_plasma(theta)
#
#    return vf
#
#def sobol_method(T):
#
#
#    x4 = np.random.rand()
#    x5 = np.random.rand()
#    x6 = np.random.rand()
#    x7 = np.random.rand()
#
#    u = -T*log(x4*x5*x6)
#    n = -T*log(x4*x5*x6*x7)
#
#    if n*n - u*u < 1:
#        return sobol_method(T)
#
#    #now we have valid u = abs(u_i)
#    x1 = np.random.rand()
#    x2 = np.random.rand()
#    #x3 = np.random.rand()
#
#    vx = u*(2*x1 -1)
#    vy = 2*u*sqrt(x1*(1-x1))
#
#    #3d treatment
#    #vy = 2*u*sqrt(x1*(1-x1))*cos(2*pi*x2)
#    #vz = 2*u*sqrt(x1*(1-x1))*sin(2*pi*x2)
#    vz = 0.0
#
#    return vx,vy,vz,u
#
##equation 6 from Swisdak
#def f1(p, A):
#    return p*p*exp(-A*p*p/(1+sqrt(1+p*p)))
#
##Derivative
##def f1p(p, A):
#
##Mode
#def f1m(A):
#    return sqrt((2.0/A/A)*(1+sqrt(1+A*A)))
#
#
#
##Rejection sampling from Swisdak 2013
#def rejection_sampling(A):
#    pm = f1m(A) #mode
#
#    #root finding
#    pg = np.linspace(0, 5, 20)
#    for p in pg:
#        print log(f1(p, A)/pm) + 1
#
#    
#    return 1.0
#
#
#
#
#
#
##cumulative distribution loading
#def drifting_maxwellian(beta, theta):
#    gamd = 1.0/sqrt(1.0-beta*beta)
#    pu = gamd*beta #drift 4-velocity
#    g1 = sqrt(1.0 + up*up)
#
#    #f(p||) 
#    fg1 = (1.0 + gamd*g1/theta)*exp(-(up-pu)**2/(g1*gamd + up*pu + 1.0)/theta)
#
#
#def boosted_maxwellian(beta, Gamma, theta):
#    
#    #For relativistic case we use Sobol method, inverse method otherwise
#    if theta > 0.1:
#        vx, vy, vz, u = sobol_method(theta)
#    else 
#        vx, vy, vz, u = inverse_method(theta)
#    
#    X8 = np.random.rand()
#    if -beta*vx > X8:
#        vx = -vx
#    else:
#        return drift_boost_maxwell(beta, Gamma, theta)
#
#    Gamma = 1.0/sqrt(1.0 - beta*beta) #XXX is this so?
#    vx = Gamma*vx + beta*sqrt(1.0 + u*u)
#
#    return vx, vy, vz, u




def thermal_plasma(theta):
    fmax = 1.0
    vmin = -5.0*theta
    vmax = 5.0*theta
    
    vf = vmin + (vmax-vmin)*np.random.rand()
    f = 0.5*(exp(-(vf*vf)/(2.0*theta*theta)))

    x = fmax*np.random.rand()

    if x > f:
        return thermal_plasma(theta)

    #now we have valid u = abs(u_i)
    x1 = np.random.rand()
    x2 = np.random.rand()
    #x3 = np.random.rand()

    vx = vf*(2*x1 -1)
    vy = 2*vf*sqrt(x1*(1-x1))

    #3d treatment
    #vy = 2*u*sqrt(x1*(1-x1))*cos(2*pi*x2)
    #vz = 2*u*sqrt(x1*(1-x1))*sin(2*pi*x2)
    vz = 0.0

    return vx,vy,vz


#relativistic maxwell-Juttner distribution
# theta is dimensionless temperature
def thermal_rel_plasma(theta):

    fmax = 1.0/kn(2,1.0/theta)
    vmin = -20.0*theta
    vmax = 20.0*theta
    vf = vmin + (vmax-vmin)*np.random.rand()
    
    f = exp(-sqrt(1+vf*vf)/theta)*vf*vf

    x = fmax*np.random.rand()

    if x > f:
        return thermal_rel_plasma(theta)

    return vf


#Sobol method to sample 4-velocities u = g*v
def sobol_method(T):
    x4 = np.random.rand()
    x5 = np.random.rand()
    x6 = np.random.rand()
    x7 = np.random.rand()

    u = -T*log(x4*x5*x6)
    n = -T*log(x4*x5*x6*x7)

    if n*n - u*u < 1:
        return sobol_method(T)

    return u



#Rejection sampling to sample relativistic momenta p = g*m*v

#equation 6 from Swisdak
def f1(p, A):
    return p*p*exp(-A*p*p/(1+sqrt(1+p*p)))

#derivative
def f1p(p, A):
    t1 = 2.0*p*exp(-A*p*p/sqrt(p*p+1)) 
    t2 = p*p*exp(-A*p*p/sqrt(p*p+1))
    t3 = A*p*p*p/(p*p + 1)**1.5
    t4 = 2*A*p/sqrt(p*p+1)

    return t1+t2*(t3-t4)
    #return 2.0/p - A*p/sqrt(p*p+1)

#log derivative
def logf1p(p, A):
    return 2.0/p - A*p/sqrt(p*p+1)

#Mode
def f1mode(A):
    return sqrt((2.0/A/A)*(1+sqrt(1+A*A)))

#root finding function 
def f1fm(p, A, fpmode):
    return log(f1(p, A)/fpmode) + 1


#Rejection sampling from Swisdak 2013
def rejection_sampling(T):
    A = 1.0/T
    #print "A=", A, "T=", 1.0/A

    pmode = f1mode(A) #mode
    fpmode = f1(pmode, A) #mode
    #print "pmode=",pmode
    #print "f(pmode)=", fpmode

    if T < 1.0: 
        pmax = 5.0/sqrt(A) #non relativistic expansion
    else:
        pmax = 12.0/A #relativistic expansion
    #print "pmax=", pmax

    pp0 = scipy.optimize.brentq(f1p, 1.0e-10, pmax, args=(A))
    #print "zero of D=", pp0
    #print "f(p0)=", f1(pp0, A)

    #root finding
    #pg = np.linspace(1.0e-10, pmax, 40)
    #for p in pg:
    #    print p, log(f1(p, A)/pmode) + 1
    #    #print p, f1p(p, A)

    pmin = 1.0e-10
    #print "start", f1(pmin, A), log(f1(pmin, A)/fpmode) +1
    #print "mid", f1(pmode, A), log(f1(pmode, A)/fpmode) +1
    #print "stop", f1(pmax, A), log(f1(pmax, A)/fpmode) +1

    pm = scipy.optimize.brentq(f1fm, pmin, pmode, args=(A, fpmode))
    pp = scipy.optimize.brentq(f1fm, pmode, pmax, args=(A, fpmode))
    #print "p- =", pm
    #print "p+ =", pp

    #now we have all the auxiliary stuff ready for the distribution
    #next lets sample with rejection method

    #lp = -f1(pp, A)/f1p(pp, A)
    #lm =  f1(pm, A)/f1p(pm, A)
    lp = -1.0/logf1p(pp, A)
    lm = 1.0/logf1p(pm, A)
    
    qm = lm/(pp-pm)
    qp = lp/(pp-pm)
    qmm =1-(qp + qm)

    X = 0.0
    while True:
        U = np.random.rand()
        V = np.random.rand()
        
        if U <= qmm:
            Y = U/qmm
            X = (1-Y)*(pm + lm) + Y*(pp - lp)
            if V <= f1(X, A)/fpmode:
                break
        elif U <= qmm + qp:
            E = -log((U-qmm)/qp)
            X = pp - lp*(1-E)
            if V <= exp(E)*f1(X, A)/fpmode:
                break
        else:
            E = -log((U-(qmm + qp))/qm)
            X = pm + lm*(1-E)
            if X > 0: #my own addition; prevents numerical underflow in log
                if V <= exp(E)*f1(X, A)/fpmode:
                    break

    #we can return X here as X = p = g*m*u, and we set m = 1 previously
    return X



#Change from u = |u| to (ux, uy, uz)
def velxy(u):

    #now we have valid u = abs(u_i)
    x1 = np.random.rand()
    x2 = np.random.rand()
    #x3 = np.random.rand()

    ux = u*(2*x1 -1)
    uy = 2*u*sqrt(x1*(1-x1))

    #3d treatment
    #uy = 2*u*sqrt(x1*(1-x1))*cos(2*pi*x2)
    #uz = 2*u*sqrt(x1*(1-x1))*sin(2*pi*x2)
    uz = 0.0

    return ux, uy, uz




#cumulative distribution loading
def drifting_maxwellian(beta, theta):
    gamd = 1.0/sqrt(1.0-beta*beta)
    pu = gamd*beta #drift 4-velocity
    g1 = sqrt(1.0 + up*up)

    #f(p||) 
    fg1 = (1.0 + gamd*g1/theta)*exp(-(up-pu)**2/(g1*gamd + up*pu + 1.0)/theta)


def boosted_maxwellian(theta, Gamma):

    #For relativistic case we use Sobol method, inverse method otherwise
    if theta > 0.2:
        #vx, vy, vz, u = sobol_method(theta)
        u = sobol_method(theta)
    else:
        #vx, vy, vz, u = rejection_sampling(theta)
        u = rejection_sampling(theta)

    #now get components
    ux, uy, uz = velxy(u)


    #We interpret this as v/c = beta
    if Gamma < 1.0:
        beta = Gamma
        Gamma = 1.0/sqrt(1.0 - beta*beta)
    else: #else as bulk lorentz factor
        beta = 1.0/sqrt(1.0 + Gamma*Gamma)
    
    #beta = 1.0/sqrt(1.0 + Gamma*Gamma)
    X8 = np.random.rand()
    if -beta*ux > X8:
        ux = -ux

    #Gamma = 1.0/sqrt(1.0 - beta*beta) #XXX is this so?
    #beta = 1.0/sqrt(1.0 + Gamma*Gamma)
    #Gamma = 1.0/sqrt(1.0-beta*beta)
    ux = Gamma*(ux + beta*sqrt(1.0 + u*u))

    u = sqrt(ux*ux + uy*uy + uz*uz)

    return ux, uy, uz, u




T = 1.0e-5
#T = 2.0
A = 1.0/T
print "T=", T, "A=", A

#beta = 1.0 - 1.0e-5
#beta = 0.99
#Gamma = 1.0/sqrt(1.0-beta*beta)
Gamma = 0.5
beta = 1.0/sqrt(1.0 + Gamma*Gamma)
print "beta=", beta, "Gamma=", Gamma

#X = drifting rejection(A, beta)
#X = rejection_sampling(A)



N = 1000
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

    ux, uy, uz, u = boosted_maxwellian(T, Gamma)
    #p = u/sqrt(1.0-u*u)
    n2[n] = ux



#plot
#print n2

#ax1.hist(n1, 100, color="black", alpha=0.3)
ax1.hist(n2, 100, color="red"  , alpha=0.3)
#ax1.hist(n3, 100, color="blue" , alpha=0.3)



fname = 'maxwells.pdf'
savefig(fname)
