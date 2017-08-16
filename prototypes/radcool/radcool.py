import numpy as np
import math
from pylab import *
import scipy
import os, sys


# Computational constants
me = 9.1093826e-28  # g
qe = 4.803250e-10   # cgs charge units
c = 2.99792458e10   # cm/s
#r_0 = 2.817940325e-13 # cm
h_pc = 6.6260693e-27        # Planck constant
sigma_T=6.65245873e-25        # Thomson cross-section



def ReBin(lnx_old,Cm_old,lnx):

    
    nn = len(lnx)                  #new photon grid defined between the 
                                   # points of the one used in program 
    nold = len(lnx_old)            #old photon grid defined in between the 
                                   # points of the precalculated one
    lnx_h=np.zeros(nn+1)
    lnx_old_h=np.zeros(nold+1)
    ind=np.zeros(nn+1,dtype=np.int)
    Cm=np.zeros(nn)

    d_lnx = (lnx[nn-1] - lnx[0])/(nn - 1)
    lnx_h = lnx - 0.5*d_lnx                  # was: lnx_h(1:nn) = lnx(1:nn) - 0.5*d_lnx
    lnx_nn=lnx[nn-1] + 0.5*d_lnx
    lnx_h=np.append(lnx_h,[lnx_nn]) # was: lnx_h(nn+1) = lnx(nn) + 0.5*d_lnx

    d_lnx_old = (lnx_old[nold-1] - lnx_old[0])/(nold - 1)
    lnx_old_h = lnx_old - 0.5*d_lnx_old
    lnx_old_nold=lnx_old[nold-1] + 0.5*d_lnx_old
    lnx_old_h=np.append(lnx_old_h,[lnx_old_nold])  # was lnx_old_h(np+1)=lnx_old(np) + 0.5*d_lnx_old

 
    if d_lnx < 2.0*d_lnx_old:
        print "ReBin: the resolution of the original grid is too small. d_lnx_old={0}, d_lnx={1}".format(d_lnx_old, d_lnx)
    
    for i in range(nn+1):                # lnx_neu
        lnxi = lnx_h[i]
        if lnxi < lnx_old_h[0]:  ind[i] = -1
        elif lnxi > lnx_old_h[nold-1]:  ind[i] = -2
        elif lnxi == lnx_old_h[0]:  ind[i] = 0
        elif lnxi == lnx_old_h[nold]:  ind[i] = -2   #is this condition ever satisfied? # ind(in) = np+1    
                                                     #If this is in, we would need Cm_old(np+1), 
                                                     # which doesn't exist.
        else:
            nmax = nold+1
            nmin = -1
            while (nmax-nmin) > 1: #locates two points of the old, 
                                   # precalculated grid, between which there is a point of the new grid 
                noldh = (nmin + nmax)//2
                if lnxi > lnx_old_h[noldh]: nmin = noldh
                elif lnxi < lnx_old_h[noldh]: nmax = noldh
                elif lnxi == lnx_old_h[noldh]:  nmin = noldh; break
            ind[i] = nmin            # assigns ind(i) minimal point in the old grid which contributes 
                                     #  to the new grid

    for i in range(1,nn+1):
        if (ind[i] < ind[i-1]) and (ind[i] != -2) and (ind[i-1] != -1):
            print "ReBin: error in determining indices. i={0}, ind(i-1)={1},ind(i)={2}".format(i,ind(i-1),ind(i))
            sys.exit() 

        Csum = 0.0            # For Cm(in-1)
        if ind[i] > -1:
            if ind[i-1] != -1:
                if (ind[i]-ind[i-1]) > 1:
                    for ii in range(ind[i-1]+1,ind[i]): # created range=ind[i-1] ... ind[i]-1
                        Csum += Cm_old[ii]*d_lnx_old

                du = lnx_h[i] - lnx_old_h[ind[i]]
                dd = lnx_old_h[ind[i-1]+1] - lnx_h[i-1]
                Csum += Cm_old[ind[i]]*du + Cm_old[ind[i-1]]*dd
            elif ind[i-1] == -1:
                if ind[i] > 0:
                    for ik in range(0,ind[i]):      ## created range=0 ... ind[i]-1 
                        Csum += Cm_old[ik]*d_lnx_old
                
                du = lnx_h[i] - lnx_old_h[ind[i]]
                Csum += Cm_old[ind[i]]*du
        elif (ind[i] == -2) and (ind[i-1] > -1):
            if ind[i-1] < nold:
                for ij in range(ind[i-1]+1, nold): # created range=ind[i-1]+1 ... nold 
                                                   # Implicitly assumes that Cm = 0 outside the old grid
                    Csum += Cm_old[ij]*d_lnx_old  

            dd = lnx_old_h[ind[i-1]+1] - lnx_h[i-1]
            Csum = Csum + Cm_old[ind[i-1]]*dd
        else:
            Csum = 0.0
        Cm[i-1] = Csum/d_lnx
    
    return Cm




# Initialisation of synchrotron emissivities on the electron/photon grid
def Synch_emis(lnx,lnz,Bfield):
    '''
    Routine reads the synchrotron emissivities from emis1.dat 
    and rebins them to the current photon grid.

    NB: the electron grid in the main program and the datafile has to match!

    1250  format(2d13.5)            # For synch from Juri
    '''
    
    # Cyclotron (electron) frequency/energy; xb = h\nu_{B}/mc^2 = eB/2\pimc * h/mc^2 
    xb = qe*Bfield/(2.0*pi*me*c)*h_pc/(me*c**2)        
    gd_const = 4.0*sigma_T/(3.0*me*c)*Bfield**2.0/(8.0*pi)

    i_m_ph = len(lnx)
    i_m = len(lnz)

    # The number of pre-calculated points in photon energy for 
    # every electron momentum (in emis1.dat)
    nx = 10000        

    d_lnx = (lnx[i_m_ph-1] - lnx[0])/(i_m_ph - 1)
    d_lnz = (lnz[i_m-1] - lnz[0])/(i_m - 1)

    pmom=np.zeros(i_m - 1)
    emigam=np.zeros((nx, i_m - 1))
    xx=np.zeros(nx)
    Cm_old=np.zeros(nx)
    lnx_old=np.zeros(nx)
    Cm=np.zeros(i_m_ph)

  
    f=open('precalc/emis1.dat','r')# A 199x10000 (z,x/xb) matrix of synch emissivities (old range)

    for k in range(len(lnz)-1):

        #reading the momentum of emitting particle
        pmom[k]=f.read(14)  #format (120,1250) pmom(k);  

        # reading photon energies and corresponding emissivities, 
        # for the specific particle momentum 
        for i in range(nx):
            xx[i] = f.read(13); emigam[i, k] = f.read(14)  
    f.close()


    # introducing ('old') energy grid, which is xb times the precalculated one
    # the precalculated values were computed for photon grid with units of xb,
    # thus one needs to multiply by the xb to come to ordinary units (?photon energy in mc^2?)
    lnx_old = log(xx*xb)  

    for k in range(len(lnz)-1):        # iteration over len(lnz)-1 elements
        z = pmom[k] # momentum of particle from the precalculated file

        # relative difference between the momentum on the grid and the precalculated momentum
        err = (z - exp(lnz[k]+0.5*d_lnz))/exp(lnz[k]+0.5*d_lnz)    
        
        # note that the momenta on the grid are supposed to be smaller by half-step,
        # i.e. z_precalc = z_grid + \delta z_grid /2. (?) This is for better cooling computations

        if err > 1e-4:
            print 'CSmatrix_j: the arrays pmom and lnz do not match. lnz(k), pmom(k), err', exp(lnz[k]+0.5*d_lnz), pmom[k], err
            sys.exit()


        # factor x is to rebin the em. ENERGY (not number) per log(x)
        Cm_old = gd_const*z*z/xb*emigam[:,k]*exp(lnx_old[:])    
        Cm = ReBin(lnx_old,Cm_old,lnx)    

        # CSmh = dN/(dlnx*dt) = 1/h_pc * dE/(dnu*dt) = dE/(d(h*nu)*dt)
        CSmh[:,k] = Cm[:]/exp(lnx[:])                

    for k in range(len(lnz)-1):
        z = pmom[k]
        gdsum = 0.0e0
        for i in range(i_m_ph):
            gdsum += CSmh[i,k]*exp(lnx[i]) # replace with matrix multiplication
        gdsum *= d_lnx

        if gdsum >= 1e-30:

            # factor gd_const*z*z/gdsum tells, how much normalization of the present synchrotron
            # emissivity coefficient differs from the theoretical one
            CSmh[:,k] *= gd_const*(z**2)/gdsum            
            if (gd_const*(z**2)/gdsum > 1.01) or (gd_const*(z**2)/gdsum < 0.99):
                if (gd_const*(z**2)/gdsum) > 1e30: 
                    print 'CSmh renorm coeff is gt 1e30:', gd_const,z,gdsum
                    sys.exit()

    return CSmh


# Calculates electron radiative cooling only, 
# averaged over pitch-angles and integrated over photon directions
# Needs electron distribution function and magnetic field
def El_synch_cool(fze,lnz,i_m,dt,Bfield):
    
    eps=1.0e-14

    # photon and electron grids (both are logarithmic)
    d_lnz = (lnz[i_m-1] - lnz[0])/(i_m - 1)
    

    d_g=np.zeros(i_m-1)
    for i in xrange(i_m-1):
        z = exp(lnz[i])
        g = (z*z + 1.0)**0.5
        zu = exp(lnz[i] + d_lnz)
        gu = (zu*zu + 1.0)**0.5
        d_g[i] = (zu - z)*(zu + z)/(gu + g)

    xb = qe*Bfield/(2.0*pi*me*c)*h_pc/(me*c**2)    # Cyclotron energy

    ################################################## 
    # Focker-Planck for cooling (equal to \dot\gamma_s)
    A_half = np.zeros(i_m-1)
    A_half = - 4.0*sigma_T*(Bfield**2/8.0/pi)*zvec**2*d_lnz/(3.0*me*c)/d_g


    ################################################## 
    # Calculating the matrix of the linear system
    M_el = np.zeros((i_m, i_m))
    for i in range(i_m):
        for i_pr in range(i-1,i+2):
            # tridiagnoal terms
            if i == 0:
                if i_pr == 0:         M_el[i,i_pr] += A_half[i]
            elif i == i_m-1:
                if i_pr == (i - 1):   M_el[i,i_pr] += - A_half[i-1]
            else:
                if i_pr == (i - 1):   M_el[i,i_pr] += - A_half[i-1]
                elif i_pr == i:       M_el[i,i_pr] += A_half[i]
                    

    # calculating matrices entering electron equation
    Mf_el=np.zeros(i_m)
    Mf_el=np.matmul(M_el,fze)

    B_el=fze/dt - c_CN*Mf_el # - c_CN*fze/t_esc 
    M_el[:,:] = (1.0 - c_CN)*M_el[:,:]
    for k in range(i_m):
        M_el[k,k] += 1.0/dt
    fz=np.linalg.tensorsolve(M_el, B_el)

    return fz





# Calculates electron cooling, radiative spectra. Used every timestep
# Needs electron distribution function, B-field
def El_evolve(fx,fze,lnx,lnz,i_m,i_m_ph,dt,Bfield,CSmh):
    
    eps=1.0e-14

    # photon and electron grids (both are logarithmic)
    d_lnx = (lnx[i_m_ph-1] - lnx[0])/(i_m_ph - 1)
    d_lnz = (lnz[i_m-1] - lnz[0])/(i_m - 1)
    d_lnz2 = d_lnz**2
    

    d_g=np.zeros(i_m-1)
    for i in xrange(i_m-1):
        z = exp(lnz[i])
        g = (z*z + 1.0)**0.5
        zu = exp(lnz[i] + d_lnz)
        gu = (zu*zu + 1.0)**0.5
        d_g[i] = (zu - z)*(zu + z)/(gu + g)

    xb = qe*Bfield/(2.0*pi*me*c)*h_pc/(me*c**2)    # Cyclotron energy

    # calculating synchrotron cooling/heating/diffusion coefficients for electron equation
    B_s_diff=np.zeros(i_m-1)
    
    for i in range(i_m-1):
        z = exp(lnz[i] + 0.5*d_lnz)
        gamma = sqrt(1.0 + z*z)

        #Locating the index of the lower boundary of integration
        lxstar = log(xb*(1.0+eps)/(gamma + z)) # argument "infinitely" close to xstar
        xstar = exp(lxstar)   # lowest energy of the non-zero emission 
                              # coefficient in first cyclotron harmonic

        if xstar <= exp(lnx[0]):
            jmin = 0
        elif xstar >= exp(lnx[i_m_ph-1]):
            B_s_diff[i] = 0.0
            continue
        else:
            for j in range(1, i_m_ph):
                if  lnx[j] >= lxstar and lnx[j-1] <= lxstar: 
                    jmin = j-1         # xmin is between the first and second integration point
        
        #integral over synchrotron emissivity times photon distribution
        Hint = 0.0
        for j in range(jmin, i_m_ph):
            x = exp(lnx[j])
            Hint = Hint + fx[j]*CSmh[j,i]*d_lnx/x
        Hfn = Hint*h_pc**2/(8.0*pi*(me*c)**3)

        # Vector of el. momenta in between gridpoints
        zvec=np.zeros(i_m-1)
        zvec[0:i_m-1] = exp(lnz[0:i_m-1] + 0.5*d_lnz)            

        # Diffusion term in the equation
        B_s_diff[i] = d_lnz*gamma*Hfn/(z**2*d_g[i])                

        # Heating term due to self-absorption
        #B1_synch_h[i] = 3.0*B1s_diff[i]                    

        if Hfn < 0.0:
            print "Electron_evol: Hfn.lt.0, Hfn = ",Hfn
            print "Electron_evol: fx = ", fx
            sys.exit()

    #gammadot_s=np.zeros(i_m-1)
    gammadot_s = 4.0*sigma_T*(Bfield**2/8.0/pi)*zvec**2*d_lnz/(3.0*me*c)/d_g     # El.cooling term, \dot{\gamma_s}

    A_half = np.zeros(i_m-1)
    B_half = np.zeros(i_m-1)
    
    # cooling and heating due to self-absorption
    ################################################## 
    #-B1_compt(1:i_m-1)+Bcc_Coul(1:i_m-1)+Bcc_adiab(1:i_m-1)+Bcc_heat(1:i_m-1) 
    # Synch cool./heat. + Compt + Coul + Adiab cooling + "external" heating (070309)
    A_half = - gammadot_s + 3.0*B_s_diff  

    # diffusion
    ################################################## 
    # + B1c_diff(1:i_m-1) + Ccc_Coul(1:i_m-1) + Ccc_heat(1:i_m-1)        
    # Sign different from Chang&cooper
    B_half =  B_s_diff                   

    
    # Determining the shifting differencing scheme according to Chang & Cooper
    d=np.zeros(i_m-1)

    for i in range(i_m-1):
        if B_half[i] == 0.0: 
            if A_half[i] > 0.0: d[i] = 1.0
            else:  d[i] = 0.0
        else:
            w = -d_lnz*A_half[i]/B_half[i]
            d[i] = 1.0/w - 1.0/(exp(w) - 1.0)    # delta in Chang&Cooper


    ################################################## 
    # Calculating the matrix of the linear system
    M_el = np.zeros((i_m, i_m))
    for i in range(i_m):
        for i_pr in range(i-1,i+2):
            # tridiagnoal temrs
            if i == 0:
                if i_pr == 0:         M_el[i,i_pr] += d[i]*A_half[i]/d_lnz + B_half[i]/d_lnz2
                elif i_pr == (i + 1): M_el[i,i_pr] += (1.0 - d[i])*A_half[i]/d_lnz - B_half[i]/d_lnz2
            elif i == i_m-1:
                if i_pr == (i - 1):   M_el[i,i_pr] += - d[i-1]*A_half[i-1]/d_lnz - B_half[i-1]/d_lnz2
                elif i_pr == i:       M_el[i,i_pr] += - (1.0 - d[i-1])*A_half[i-1]/d_lnz + B_half[i-1]/d_lnz2
            else:
                if i_pr == (i - 1):   M_el[i,i_pr] += - d[i-1]*A_half[i-1]/d_lnz - B_half[i-1]/d_lnz2
                elif i_pr == i:       M_el[i,i_pr] += (d[i]*A_half[i] - (1.0 - d[i-1])*A_half[i-1])/d_lnz \
                                                       + (B_half[i] + B_half[i-1])/d_lnz2
                elif i_pr == (i + 1): M_el[i,i_pr] += (1.0 - d[i])*A_half[i]/d_lnz - B_half[i]/d_lnz2
                    

    # calculating matrices entering electron equation
    Mf_el=np.zeros(i_m)
    Mf_el=np.matmul(M_el,fze)

    B_el=fze/dt - c_CN*Mf_el # - c_CN*fze/t_esc 
    M_el[:,:] = (1.0 - c_CN)*M_el[:,:]
    for k in range(i_m):
        M_el[k,k] += 1.0/dt
    fz=np.linalg.tensorsolve(M_el, B_el)
    

    # check if evolution is going too fast

    # Test of total energy gain
    # En_gain_el = 0.0
    # for k in range(i_m):
    #     z = exp(lnz[k])
    #     gamma = (z*z + 1.0)**0.5
    # unupdated d_t has to be used. Wrong if new d_t isn't d_t*gt
    #     En_gain_el += (fz[k] - fz_prev[k])*gamma*d_lnz/dt		
    # print 'Electron energy gain', En_gain_el

    Em_sum = 0.0
    Abs_sum = 0.0
    for k in range(i_m-1):
        t1 = d_g[k]*3.0*B_s_diff[k]/d_lnz
        t2 = d_g[k]*B_s_diff[k]/d_lnz
        tmp1=(1.0-c_CN)*fz[k+1] + c_CN*fze[k+1]
        tmp2=(1.0-c_CN)*fz[k] + c_CN*fze[k]
        Abs_sum += d_lnz*t1*( (1.0 - d[k])*tmp1 + d[k]*tmp2) - t2*(tmp1 - tmp2)
        Em_sum += d_g[k]*gammadot_s[k]*((1.0 - d[k])*tmp1 + d[k]*tmp2)
    Engain_el_Synch = Abs_sum - Em_sum

    return [fz, d]
 



# Calculates photons distribution. Used every timestep
# Needs particle distribution function, B-field
def Ph_evolve(fx,fz,lnx,lnz,i_m,i_m_ph,dt,Bfield,CSmh,d):

    xb = qe*Bfield/(2.0*pi*me*c)*h_pc/(me*c**2) # Cyclotron energy

    eps=1.0e-8

    # photon and electron grids (both are logarithmic)
    d_lnx = (lnx[i_m_ph-1] - lnx[0])/(i_m_ph - 1)
    d_lnz = (lnz[i_m-1] - lnz[0])/(i_m - 1)
    d_lnz2 = d_lnz**2
    zmin = exp(lnz[0])
    zmax = exp(lnz[i_m-1])

    d_g = np.zeros(i_m-1)
    for i in xrange(i_m-1):
        z = exp(lnz[i])
        g = (z*z + 1.0)**0.5
        zu = exp(lnz[i] + d_lnz)
        gu = (zu*zu + 1.0)**0.5
        d_g[i] = (zu - z)*(zu + z)/(gu + g)

    M_ph = np.zeros((i_m_ph, i_m_ph))
    B_ph = np.zeros(i_m_ph)

    # Determining the absorption and emission terms in the equation
    for i in range(i_m_ph):  
        x = exp(lnx[i])
        lzstar = log(0.5*(1.0+eps)*(xb/x - x/xb))  # argument "infinitely" close to zstar
        zstar = exp(lzstar)

        jsum = 0.0
        ksum = 0.0
        for j in range(i_m-1):
            z = exp(lnz[j] + 0.5*d_lnz)
            gamma = sqrt(z*z + 1.0)
            gPz2 = CSmh[i,j]*gamma/(z**2)
            fz12 = (1.0 - d[j])*fz[j+1] + d[j]*fz[j]
            dfz = (fz[j+1] - fz[j])/d_lnz
            jsum = jsum + fz12*CSmh[i,j]
            ksum = ksum + (3.0*gPz2*fz12 - gPz2*dfz)*d_lnz

        #emission term
        B_ph[i] = jsum*d_lnz/h_pc            

        #absorption term
        M_ph[i,i] = ksum/(x**2)
        if zstar >= zmax:
            B_ph[i] = 0.0
            M_ph[i,i] = 0.0

    #-------------------------------------------------- 
    # test emissivity
    #Bcr=2.0*pi*(me**2)*(c**3)/qe/h_pc
    #b=Bfield/Bcr
    #Ub=Bfield**2/(8.0*pi)
    #
    #for k in range(i_m_ph):
    #    if CSmh[k,111]>0: print exp(lnx[k]-log(b)), CSmh[k,111]*b/(4.0*sigma_T*Ub/(3.0*me*c))/h_pc
    #-------------------------------------------------- 


    # (self) Absorption term, c*k
    # The power of c is 3 instead of 4 because in the equation we have c*k
    M_ph = M_ph*h_pc**2.0/(8.0*pi*(me*c)**3.0)


    # Calculating the matrix of the linear system 
    # M_ph = M_ph_s #+ M_ph_c + M_ph_pp    # Synchrotron coefficients into matrix
    Mf_ph=np.zeros(i_m_ph)
    for i in range(i_m_ph):
      Mf_ph[:] += M_ph[:,i]*fx[i]
    B_ph = B_ph + fx/dt - c_CN*Mf_ph #- c_CN*fx/t_esc
    for k in range(i_m_ph):
        M_ph[k,k] = M_ph[k,k]*(1.0 - c_CN) + 1.0/dt  #+ (1.0 - c_CN)/t_esc
    fx_new=linalg.tensorsolve(M_ph, B_ph)


    #clip
    for i in range(i_m_ph):
        if fx[i] < 0.0: fx[i]=0.0
    
    #energy gains
    Engain_ph = 0.0
    for i in range(i_m_ph):
        x = exp(lnx[i])
        for  j in range(i_m_ph):
            tmp = (1.0 - c_CN)*fx_new[j] + c_CN*fx[j]
            Engain_ph = Engain_ph + M_ph[i,j]*tmp*x
    Engain_ph = - Engain_ph*d_lnx
 
    Em_sum = 0.0
    for k in range(i_m_ph):
        x = exp(lnx[k])
        Em_sum += x*B_ph[k]*d_lnx

    Engain_ph_Synch = Engain_ph + Em_sum

    return fx_new


def plot_pdf_el(ax, lnz, fz, fz_init, xlabel='', color='blue'):
    
    ax.cla()
    ax.set_xlabel(xlabel)
    ax.set_xlim(3e-2, 1e3)
    ax.set_ylim(1e-4, 3e0)
    ax.set_yscale('log')
    ax.set_xscale('log')

    #ax.set_xlim(0.0, 20.0)


    ax.plot(exp(lnz), fz, color=color, alpha=0.8)
    ax.plot(exp(lnz), fz_init, color='red', alpha=0.8)

    return ax


def plot_pdf_ph(ax, lnx, fx, fx_init, xlabel='', color='blue'):
    
    ax.cla()
    ax.set_xlabel(xlabel)
    ax.set_ylim(1e6, 5e35)
    ax.set_yscale('log')
    ax.set_xscale('log')


    #ax.set_ylim(1e30, 1e35)
    #ax.set_xlim(0.0, 1.0)

    ax.plot(exp(lnx), fx, linestyle='solid', marker='.',  color=color, alpha=0.8)
    ax.plot(exp(lnx), fx_init, "r-", alpha=0.8)

    return ax





##################################################
##################################################
##################################################


i_m_ph=200                  # Number of photon grid points
x_min=1.0e-11                # Minimal photon energy, x = h\nu / mc^2
x_max=9772.3722e0            # Maximal photon energy
i_m=200                      # Number of electron/positron grid points 
z_min = 1.0e-4               # Minimal electron/positron momentum, z = p/mc
z_max = 1.0e4                # Maximal electron/positron momentum
dt=1.0e-4                    # Time step

R=1e7       # Size of the medium
t_esc=R/c   # Escape (light-crossing) time
tau=1.0     # initial Thomson optical depth
Bfield=1e5                   # Magnetic field, [G]


lnx=np.zeros(i_m_ph)
lnz=np.zeros(i_m)
lnx_min = np.log(x_min)
lnx_max = np.log(x_max)


# Array of photon energy logarithms
for j in range(i_m_ph):
  lnx[j] = lnx_min + j*(lnx_max - lnx_min)/(i_m_ph - 1.0)  
d_lnx = (lnx[i_m_ph-1] - lnx[0])/(i_m_ph - 1.0)

# Array of electron momentum logarithms
lnz_min = np.log(z_min)
lnz_max = np.log(z_max)
for j in range(i_m):
  lnz[j] = lnz_min + j*(lnz_max - lnz_min)/(i_m - 1.0)        
d_lnz = (lnz[i_m-1] - lnz[0])/(i_m - 1.0)


# Reading synchrotron emissivities
CSmh=np.zeros((i_m_ph, i_m-1))
CSmh=Synch_emis(lnx,lnz,Bfield)  # Reading the precalculated emissivities and rebinning
CSmh = CSmh*h_pc                 # Transforming from dN_ph/dlnx representation to dE_ph/d\nu representation.


fze=np.zeros(i_m)
fz=np.zeros(i_m)
fx=np.zeros(i_m_ph)
d=np.zeros(i_m-1)


##################################################
# Path to be created
path = "out"
if not os.path.exists(path):
    os.makedirs(path)


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


z=exp(lnz)
gamma=sqrt(z**2 + 1.0)
#fze=gamma**(-4.0)*exp(lnz)**2  # df/d\gamma = gamma^-3; df/dlnz = df/d\gamma * z^2/gamma
fze=gamma**(-3.0)*exp(lnz)**2  # df/d\gamma = gamma^-3; df/dlnz = df/d\gamma * z^2/gamma
#fze = z**3.0*exp(-z**2/(0.1*(gamma+1.0))) # initial electron distribution

fz_int=np.trapz(fze, dx=d_lnz)

# Setting init. Th. opt. thickness eq. to the equil. value. 
# Also needed to account for background el. for comparing with Coppi
fze = fze*tau/(fz_int*sigma_T*R)		

# initial photon distribution
#fx[:]=0.0            
x=exp(lnx)
fx=8.0*pi*(me*c*x)**3.0/(h_pc**3.0*(exp(x/2.0e-5) - 1.0))




nsteps=100

#starting distributions
fz_init=fze
Lx_init=fx*exp(lnx)*4*pi*R**3/(3*t_esc)/1.22e6


for step in range(nsteps+1):
    print "Timestep = ", step+1
    
    # Crank-Nicolson coefficient, if set to vary, then i is the timestep
    c_CN=0.5*(1.0 - exp(-1.0*(step+1)**4/1.0e5))   


    # Evolving electron distribution
    fze_new, d = El_evolve(fx,fze,lnx,lnz,i_m,i_m_ph,dt,Bfield,CSmh)        

    # Evolving positron/ion distribution
    # fzp, d = El_evolve(fx,fzp,lnx,lnz,i_m,i_m_ph,dt,Bfield)    

    print c_CN
    print d

    #fz=fze # + fzp    # Sum up over all emitting particle species


    fz_int=np.trapz(fze_new, dx=d_lnz)

#    print "El distribution function integral", fz_int*sigma_T*R
    
#     En_gain_el = 0.0
#     for k in range(i_m):
#         z = exp(lnz[k])
#         gamma = (z*z + 1.0)**0.5

#         unupdated d_t has to be used. Wrong if new d_t isn't d_t*gt
#         En_gain_el += (fz[k] - fz_prev[k])*gamma*d_lnz/dt		
#     print 'Electron energy gain per sec, in mc^2 units, per dV', En_gain_el
#    fx_prev=fx


    fx_new = Ph_evolve(fx,fze_new,lnx,lnz,i_m,i_m_ph,dt,Bfield,CSmh,d)

#     En_gain_ph = 0.0
#     for j in range(i_m_ph):
#         x = exp(lnx[j])
#         # unupdated d_t has to be used. Wrong if new d_t isn't d_t*gt
#         En_gain_ph += (fx[j] - fx_prev[j])*x*d_lnx/dt		
#     print 'Photon energy gain per sec, in mc^2 units, per dV', En_gain_ph



    fze=fze_new
    fx=fx_new

    if (step % 20) == 0: 
        Lx=fx_new*exp(lnx)*4*pi*R**3/(3*t_esc)/1.22e6
        #print "Ph Luminosity", np.trapz(fze_new*exp(lnx)*4*pi*R**3/(3*t_esc)/1.22e6, dx=d_lnx)

        ax1 = plot_pdf_ph(ax1, lnx, Lx, Lx_init, r'$x=h\nu/m_e c^2$')
        ax2 = plot_pdf_el(ax2, lnz, fze_new*sigma_T*R, fz_init*sigma_T*R, '$z=p/mc$')

        fname = path+'/rad_'+str(step)+'.png'
        savefig(fname)

