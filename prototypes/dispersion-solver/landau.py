import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy.special import erf
from scipy.optimize import minimize_scalar
from math import isnan
from math import isinf


from dispsol import Jpole8, Jpole12
from dispsol import ES1d




plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=9)

fig = plt.figure(figsize=(3.54, 4.0)) #single column fig
#fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(2, 1, hspace=0.15)


axs = []
axs.append( plt.subplot(gs[0,0]) )
axs.append( plt.subplot(gs[1,0]) )

for ax in axs:
    ax.minorticks_on()
    ax.set_xlabel(r'$k \lambda_D$')
    #ax.set_xlim((0.0, 100.0))

axs[0].set_ylabel(r'Re{ $\omega/\omega_p$ }')
axs[1].set_ylabel(r'Im{ $\omega/\omega_p$ }')

axs[0].set_ylim((0.9,   1.5))
#axs[1].set_ylim((0.02, -0.16))
axs[1].set_ylim((-0.16, 0.02))

#logscale growth rates
#axs[1].set_ylim((1.0e-6, 1.0e0))
#axs[1].set_yscale('log')


#Langmuir wave Landau damping
qs  = np.array([1.0])
ms  = np.array([1.0])
ns0 = np.array([1.0])
Ts  = np.array([1.0])
vs0 = np.array([0.0])


wps  = np.sqrt(ns0 * qs**2.0 / ms)   #plasma frequency
vts  = np.sqrt(2.0*Ts/ms)            #thermal speed
lDeb = np.sqrt(Ts / (ns0 * qs**2.0)) #Debye length
kDeb = 1.0/lDeb                      #Debye length


print("vth= ", vts)
print("ldeB=", lDeb)
print("kDeb=", kDeb)
print("wps= ", wps)
print("vs0= ", vs0)

lDebTot = np.sqrt(1.0/np.sum(ns0/Ts))
print("lDebtot:", lDebTot)

print("vs/vth:", vs0/vts)


##################################################
# testing J-Pole expansion
#bj, cj = Jpole8()
bj, cj = Jpole12()

J = len(bj)  #number of pole expansion
S = len(ns0) #number of species

print("sum bj     : ", np.sum(bj))
print("sum bj*cj  : ", np.sum(bj*cj))
print("sum bj*cj^2: ", np.sum(bj*cj**2.0))



# visualize Langmuir wave dispersion relation
params = {'vts': vts, 
          'vs0': vs0, 
          'lDeb': lDeb, 
          'S':S, 
          'J':J}

karr = np.linspace(0.01, 0.5, 100)
warr = np.zeros((len(karr), S*J), np.complex)

for i,k in enumerate(karr):
    w = ES1d(k, params)
    warr[i,:] = w[:]

ms = 1.0
Nsol = 1
for nsol in range(Nsol):
    axs[0].plot(karr,  np.abs(np.real(warr[:,nsol])), 'k-', markersize=ms)

#axs[1].plot(karr, np.zeros(len(karr)), "r--")
for nsol in range(Nsol):
    axs[1].plot(karr, np.imag(warr[:,nsol]), 'k-', markersize=ms)


# for printing
karrp = [0.0496729413289805, #mode1
         0.099345882657961,  #2
         0.1490188239869415, #3
         0.198691765315922,  #4
         0.2483647066449025, #5
         0.298037647973883,  #6
         0.3477105893028635, #7
         0.397383530631844,  #8
         0.4470564719608245, #9
         0.496729413289805,  #10
         ]

for i,k in enumerate(karrp):
    w = ES1d(k, params)

    print("mode=",i+1)
    print("khat=",k)
    print("omeg=",w[1])
    print()



plt.subplots_adjust(left=0.18, bottom=0.09, right=0.98, top=0.95, wspace=0.0, hspace=0.0)
plt.savefig('landau.pdf')
