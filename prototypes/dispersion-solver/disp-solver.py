import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy.special import erf
from scipy.optimize import minimize_scalar
from math import isnan
from math import isinf


plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)

fig = plt.figure(figsize=(3.54, 5.0)) #single column fig
#fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(3, 1, hspace=0.15)


axs = []
axs.append( plt.subplot(gs[0,0]) )
axs.append( plt.subplot(gs[1,0]) )
axs.append( plt.subplot(gs[2,0]) )

for ax in axs:
    ax.minorticks_on()
    ax.set_xlabel(r'$k \lambda_D$')
    #ax.set_xlim((0.0, 100.0))

axs[0].set_ylabel(r'Re{ $\omega$ }')
axs[1].set_ylabel(r'Im{ $\omega$ }')



#Langmuir wave Landau damping
#qs  = np.array([1.0])
#ms  = np.array([1.0])
#ns0 = np.array([1.0])
#Ts  = np.array([1.0])
#vs0 = np.array([0.0])

#bump-on-tail
#qs  = np.array([-1.0,-1.0])
#ms  = np.array([1.0, 1.0])
#ns0 = np.array([0.9, 0.1])
#Ts  = np.array([1.0, 1.0])
#vs0 = np.array([0.0, 7.071])

#two-stream
qs  = np.array([-1.0,  -1.0])
ms  = np.array([ 1.0,   1.0])
ns0 = np.array([ 0.5,   0.5])
Ts  = np.array([ 0.5,   0.5])
vs0 = np.array([ 0.0,   13.0])
 
#Ts[1]  = 0.5*vs0[1]*vs0[1]*(1.0/5.0)**2.0 # in units of vb

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


def Jpole12():
    """ 12th order J-pole expansion of plasma zeta function """

    N = 12
    bj= np.zeros((N), np.complex)
    cj= np.zeros((N), np.complex)

    bj[0]= -0.00454786121654587 - 0.000621096230229454j;
    bj[1]=  0.215155729087593      + 0.201505401672306j;
    bj[2]=  0.439545042119629 +       4.16108468348292j;
    bj[3]=  -20.2169673323552 -       12.8855035482440j;
    bj[4]=  67.0814882450356 +        20.8463458499504j;
    bj[5]=  -48.0146738250076 +       107.275614092570j;
                                                        
    cj[0]=  -2.97842916245164 -       2.04969666644050j;
    cj[1]=  2.25678378396682 -        2.20861841189542j;
    cj[2]=  -1.67379985617161 -       2.32408519416336j;
    cj[3]=  -1.15903203380422 -       2.40673940954718j;
    cj[4]=  0.682287636603418 -       2.46036501461004j;
    cj[5]=  -0.225365375071350 -      2.48677941704753j;

    bj[6:] = np.conj(bj[0:6]);
    cj[6:] =-np.conj(cj[0:6]);

    return bj, cj

def Jpole8():
    """ 8th order J-pole expansion of plasma zeta function
        (single prec) Ronnmark 1982 """
 
    N = 8
    bj= np.zeros((N), np.complex)
    cj= np.zeros((N), np.complex)

    bj[0]=-1.734012457471826e-2 -4.630639291680322e-2j
    bj[1]=-7.399169923225014e-1 +8.395179978099844e-1j
    bj[2]=5.840628642184073     +9.536009057643667e-1j
    bj[3]=-5.583371525286853    -1.120854319126599e1j

    cj[0]=2.237687789201900-1.625940856173727j
    cj[1]=1.465234126106004-1.789620129162444j
    cj[2]=.8392539817232638-1.891995045765206j
    cj[3]=.2739362226285564-1.941786875844713j

    #double prec
    #bj[0]=  -0.0173401116032742 - 0.0463064419344598j
    #bj[1]=  -0.739917851897683 + 0.839518298070637j
    #bj[2]=  5.84063227513760 + 0.953602843950785j
    #bj[3]=  -5.58337431170864 - 11.2085508179677j
    #     
    #cj[0]=   2.23768772215616 - 1.62594103256666j
    #cj[1]=   1.46523409042510 - 1.78962030806222j
    #cj[2]=   0.839253965702731 - 1.89199521968963j
    #cj[3]=   0.273936217871668 - 1.94178704551807j

    bj[4:]= np.conj(bj[0:4]);
    cj[4:]=-np.conj(cj[0:4]);

    return bj, cj



# get J-pole expansion and check validity
bj, cj = Jpole12()
#bj, cj = Jpole8()

J = len(bj)  #number of pole expansion
S = len(ns0) #number of species

print("sum bj     : ", np.sum(bj))
print("sum bj*cj  : ", np.sum(bj*cj))
print("sum bj*cj^2: ", np.sum(bj*cj**2.0))



def ES1d(k, params):
    """ solve omega from plasma dispersion using matrix decomposition """

    if k <= 0.0:
        return np.zeros((2))
    if k > 100.0:
        return np.zeros((2))

    global bj, cj
    global J,S

    vts  = params['vts']
    vs0  = params['vs0']
    lDeb = params['lDeb']


    M = np.zeros((S*J, S*J), np.complex)

    bsj = np.zeros((S*J), np.complex)
    csj = np.zeros((S*J), np.complex)
    
    #ES1D dense matrix
    sj = 0
    for s in range(S):
        for j in range(J):
            csj[sj] = k*( cj[j]*vts[s] + vs0[s] );
            #bsj[sj] = bj[j]*cj[j]*vts[s]/( k*kDeb[s]**2.0 );
            bsj[sj] = bj[j]*cj[j]*vts[s]/(k*lDeb[s]**2.0);
            sj += 1

    for sj in range(S*J):
        M[sj,: ] = -bsj[sj];
        M[sj,sj] = -bsj[sj] + csj[sj];

    omega = np.linalg.eigvals(M)

    indx = np.argsort(-np.imag(omega) ) #sort in descending order in place
    #omega[indx] #sort in descending order in place

    #for sj in range(S*J):
    #    M[sj,: ] = -sj
    #    M[sj,sj] = sj

    #print(M)

    return omega[indx]


#return the fastest growing mode's growth rate
def fastest_growing_mode(k, params):
    omega = ES1d(k, params) 

    omega1 = np.imag(omega[0])
    if isnan(omega1):
        omega1 = 0.0 + 0.0j
    if isinf(omega1):
        omega1 = 0.0 + 0.0j

    return -omega1


def find_max_growth(params):

    #initial_guess = [x0]
    res = minimize_scalar(
            fastest_growing_mode,
            args=params,
            #bounds=(0.01, 0.4),
            #method='bounded',
            bracket=(1.0e-3, 0.4),
            method='brent',
            tol=1.0e-4
            )

    return res


params = {'vts': vts, 'vs0': vs0, 'lDeb': lDeb}

kmode = 0.2
omega = ES1d(kmode, params)
print("test solution:", omega[:3])


##################################################
#sys.exit()
# visualize Langmuir wave dispersion relation

#karr = np.linspace(0.01*lDeb[0], 1.0*lDeb[0], 100)
karr = np.linspace(0.01, 1.0, 100)
#karr = np.logspace(-3, np.log10(0.3), 100)
warr = np.zeros((len(karr), 16), np.complex)

for i,k in enumerate(karr):
    #print(k*kDeb)
    w = ES1d(k, params)
    warr[i,:] = w[:16]

ms = 1.0
Nsol = 4
for nsol in range(Nsol):
    axs[0].plot(karr,  np.real(warr[:,nsol]), '.', markersize=ms)

axs[1].plot(karr, np.zeros(len(karr)), "k--")
for nsol in range(Nsol):
    axs[1].plot(karr,  np.imag(warr[:,nsol]), '.', markersize=ms)


res = find_max_growth(params)
print(res.x)
print(-res.fun)

sys.exit()




##################################################

vbeam  = np.linspace(1.0, 15.0, 90)
grates = np.zeros((len(vbeam)))
kmodes = np.zeros((len(vbeam)))

Gmmax  = np.zeros((len(vbeam)))

for i, x in enumerate(vbeam):

    #twostream
    #qs  = np.array([-1.0,  -1.0])
    #ms  = np.array([ 1.0,   1.0])
    #ns0 = np.array([ 0.5,   0.5])
    #Ts  = np.array([ 0.5,   0.5])
    #vs0 = np.array([ 0.0,   np.sqrt(2)*x])
     
    qs  = np.array([-1.0,  -1.0,  1.0,   1.0])
    ms  = np.array([ 1.0,   1.0,  1.0,   1.0])
    ns0 = np.array([ 1.0,   1.0,  1.0,   1.0 ])
    Ts  = np.array([ 0.5,   0.5,  0.5,   0.5])
    vs0 = np.array([ 0.0,   x, 0.0, x])
    
    #vs0 *= np.sqrt(2.0) #normalize

    #bump-on-tail setup
    #alpha = 1.0e-1
    #ns0[0] = 1.0-alpha
    #ns0[1] = alpha
    #ns0[2] = 1.0-alpha
    #ns0[3] = alpha
    #ns0 /= 2.0
    ns0 /= 4.0

    wps  = np.sqrt(ns0 * qs**2.0 / ms)   #plasma frequency
    vts  = np.sqrt(2.0*Ts/ms)            #thermal speed

    lDeb = np.sqrt(Ts / (ns0 * qs**2.0)) #Debye length
    kDeb = 1.0/lDeb                      #Debye length
    
    params = {'vts': vts, 'vs0': vs0, 'lDeb': lDeb}
    res = find_max_growth(params)

    grates[i] = res.fun
    kmodes[i] = res.x

    #theoretical prediction
    alpha = ns0[1]/ns0[0]
    vb = vs0[1]
    Vb = np.sqrt(Ts[1])
    Gmmax[i] = np.sqrt(np.pi/(2.0*np.exp(1.0))) * alpha * (vb/Vb)**2.0


axs[2].set_xlabel(r'$v_b/\sqrt{\theta}$')
axs[2].set_ylabel(r'$\omega$')

axs[2].plot(vbeam, -grates, 'r-')
axs[2].plot(vbeam,  kmodes, 'b-')

axs[2].set_ylim((0.0, 0.4))

#analytics for two-stream

Gm = 0.5/np.sqrt(2) #predicted growth rate saturation
print("predicted growth rate saturation:", Gm)
axs[2].plot(vbeam, Gm*np.ones(len(vbeam)), 'k--')

#axs[2].plot(vbeam, Gmmax, 'r--')
#print(Gm/Gmmax[-1])



plt.subplots_adjust(left=0.18, bottom=0.07, right=0.98, top=0.95, wspace=0.0, hspace=0.0)
plt.savefig('dispersion.pdf')
