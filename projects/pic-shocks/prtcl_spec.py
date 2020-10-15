import numpy as np
import matplotlib.pyplot as plt

#import pymc3 as mcmc
#import emcee
from scipy.optimize import minimize
import h5py as h5

import sys, os
from parser import parse_input
import argparse

from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import colorbar

from particle_spectra import default_values
from particle_spectra import default_turbulence_values
from particle_spectra import create_spectra

#from comp_spectra import find_min
#from comp_spectra import find_max



def model_f2(theta, x):
    psi, zeta, K = theta    

    gam0 = 200.0
    gamth = 1.1
    C = 1.0
    A = C/gamth
    gamcool = ( ((3.0-psi)/A)*zeta*gam0**(2.0-psi) )**(1./(3.0-psi))

    return K*np.exp(-(x/gamcool)**(3.0-psi))

def find_min(xx, th):
    i1 = 0
    for i,x in enumerate(xx):
        if x >= th:
            i1 = i
            break
    return i1

#def model_f(theta, x):
#    psi, gamcool, K = theta
#    return K*np.exp(-(x/gamcool)**(3.0-psi))

def model_f(theta, x):
    psi, gamcool, K, p = theta
    #plaw = K*x**(-p)
    #plaw = K*x**(0)
    #plaw = K/x
    plaw = K
    return plaw*np.exp(-(x/gamcool)**(3.0-psi))


def log_likelihood(theta, x, y):
    model = model_f(theta, x)

    nw = len(x)
    sigma2 = np.ones(nw)
    #sigma2 = y
    return -0.5*np.sum((y-model)**2/sigma2 + np.log(sigma2))


def log_prior(theta):
    #psi, gamcool, K = theta    
    psi, gamcool, K, p = theta    
    if 0.0 < psi < 3.0 and 5.0 < gamcool < 40.0 and 1.0e2 < K < 1.0e5 and -2.0 < p < 4.0:
        return 0.0
    return -np.inf


def log_probability(theta, x, y, gaminj):
    lp = log_prior(theta)
    #i1 = find_min(x, gaminj)
    i1 = find_min(x, 0.8*theta[1])

    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x[i1:], y[i1:])




def fit_gamma_cool(bins, hist, color='k'):

    axs[0].plot(bins, hist, color=color, linewidth=0.7)

    return 

    #--------------------------------------------------
    # fit
    if False:
        #nll = lambda *args: -log_likelihood(*args)
        nll = lambda *args: -log_probability(*args)
        #initial = np.array([1.8, 20.0, 1.0e2, -0.8]) 
        #initial = np.array([1.5, 20.0, 7.0e3, 0.5]) 
        initial = np.array([1.7, 10.0, 2.0e3, 0.0]) 

        if False:         
            probs = np.zeros(50)
            gaminjs = np.linspace(1.0, 40.0, 50)
            for i,gaminj in enumerate(gaminjs):
                soln = minimize(nll, initial, args=(bins, hist, gaminj))
                #m_psi, m_gamcool, m_K, m_gaminj = soln.x

                probs[i] = soln.fun
                print(gaminj, np.log10(soln.fun))

            # get best gaminj
            i = np.argmax(probs)
            gaminj = gaminjs[i]
        else:
            gaminj = 10.0

        soln = minimize(nll, initial, args=(bins, hist, gaminj), method='Nelder-Mead')
        m_psi, m_gamcool, m_K, m_p = soln.x

        print("Maximum likelihood estimates: {0:.5f}".format(soln.fun))
        print("psi  = {0:.3f}".format(m_psi))
        print("gamc = {0:.3f}".format(m_gamcool))
        print("K    = {0:.3f}".format(m_K))
        print("p    = {0:.3f}".format(m_p))
        print("gami = {0:.3f}".format(gaminj))

        #i1 = find_min(bins, gaminj)
        i1 = find_min(bins, 0.8*m_gamcool)
        xx = bins[i1:]
        yy = hist[i1:]
        axs[0].plot(xx, yy, color='blue', alpha=0.6)

        yy1 = model_f(soln.x, xx)
        yy2 = model_f(initial, xx)
        axs[0].plot(xx, yy1, "r--")
        axs[0].plot(xx, yy2, "g:")

        #yy3 = m_K*xx**(-m_p)
        #axs[0].plot(xx, yy3, color='darkorange', linestyle="dotted")

def get_top_particles(
        ntop,
        info, 
        isp  = None,
        read_file = False,
        ):

    # test particle spectra
    if True:
        fname = info['particle_file'] + "_" + str(info['lap']) + ".h5"
        f5F = h5.File(fname,'r')
        nx = f5F['Nx'].value
        ny = f5F['Ny'].value
        nz = f5F['Nz'].value
        #print("reshaping 1D array into multiD with ({} {} {}) {}".format(nx,ny,nz, info['lap']))

        ux   = np.reshape( f5F['vx'][:]   ,(nx,-1)) 
        uy   = np.reshape( f5F['vy'][:]   ,(nx,-1)) 
        uz   = np.reshape( f5F['vz'][:]   ,(nx,-1)) 
        #wgt  = np.reshape( f5F['wgt'][:]  ,(nx,-1)) 
        gam = np.sqrt(1.0 + ux*ux + uy*uy + uz*uz).flatten()
        
    
    gam = np.sort(gam)
    #gam = gam[-10:] #drop 10 most energetic ones out
    #return gam[-ntop:]
    return gam[-ntop:]

def top_particles_id(
        ntop,
        info,
        topgam,
        isp = None,
        read_file = False,
        ):
    
    if True:
        fname = info['particle_file'] + "_" + str(info['lap']) + ".h5"
        f5F = h5.File(fname,'r')
        nx = f5F['Nx'].value
        ny = f5F['Ny'].value
        nz = f5F['Nz'].value
        i_d = f5F['id'].value
        proc = f5F['proc'].value
        
        ux   = np.reshape( f5F['vx'][:]   ,(nx,-1)) 
        uy   = np.reshape( f5F['vy'][:]   ,(nx,-1)) 
        uz   = np.reshape( f5F['vz'][:]   ,(nx,-1))
        gam = np.sqrt(1.0 + ux*ux + uy*uy + uz*uz).flatten()
        
        top_10_array_id = []
        top_10_array_proc = []
        print(topgam)
        print(gam)
        check = 0
        topcheck = 0
        for check in range(len(gam)):
            for topcheck in range(ntop):
                if gam[check] == topgam[topcheck]:
                    top_10_array_id.append(i_d[check])
                    top_10_array_proc.append(proc[check])
    print("The id and proc of top {} particles is:".format(ntop))
    print(top_10_array_id)
    print(top_10_array_proc)
    return top_10_array_id, top_10_array_proc

def lap2time(lap, conf):
    L = conf.Ny*conf.NyMesh
    l0 = L
    t0 = l0/conf.cfl

    print("lap: {} t: {} t0: {} L: {} l_0: {}".format(lap,lap/t0, t0, L, l0))
    return lap/t0


if __name__ == "__main__":

    fig = plt.figure(1, figsize=(3.487, 2.5))

    plt.rc('font',  family='sans-serif')
    #plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.3)
    #gs.update(wspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    
    for ax in axs:
        ax.minorticks_on()

    axs[0].set_yscale('log')
    axs[0].set_xscale('log')

    axs[0].set_ylim((1.0e0, 3e4))
    axs[0].set_xlim((1, 100.0))
                           
    axs[0].set_xlabel(r"$\gamma$")
    axs[0].set_ylabel(r"$dN/d\ln(\gamma)$")

    #--------------------------------------------------
    # command line driven version
    conf, fdir, args = parse_input()
    fdir += '/'
    print("processing {}".format(fdir))
    #fname_P = "particles"
    fname_P = "test-prtcls"


    # if lap is defined, only process that one individual round
    if not(args.lap == None):
        info = {}
        info['lap'] = args.lap
        info['fdir'] = fdir
        #info['particle_file'] = fdir+'restart/'+fname_P
        info['particle_file'] = fdir+fname_P
        bins, hist = create_spectra('gamma', info, conf, read_file=False, no_plots=True)

        #skim of small bin counts
        #fi = np.where(hist > 20)
        #bins = bins[fi]
        #hist = hist[fi]
        fit_gamma_cool(bins, hist)


    # else process every file that there is
    else:
        norm = Normalize(vmin=0.0, vmax=10.)
        cmap = get_cmap('Spectral')

        files_P = []
        #for lap in range(conf.restart, conf.Nt+1, conf.restart):
        #for lap in range(conf.interval, conf.Nt+1, conf.interval):

        merge_step = 1
        for lap in range(0, conf.Nt+1, merge_step*conf.interval):
            #files_P.append( fdir+'restart/'+fname_P)
            files_P.append( fdir+'restart/'+fname_P)
    
            #info = build_info(fdir, lap)
            info = {}
            info['lap'] = lap
            info['fdir'] = fdir
            #info['particle_file'] = fdir+'restart/'+fname_P
            info['particle_file'] = fdir+fname_P

            if not(os.path.isfile(info['particle_file'] + "_" + str(info['lap']) + ".h5")):
                continue

            bins, hist = create_spectra('gamma', info, conf, read_file=False, no_plots=True)

            laptime = lap2time(lap,conf)
            topgams = get_top_particles(10, info)
            print("top10 t={}: {}, {}".format(laptime, np.mean(topgams), np.std(topgams)))
            


            #if laptime > 3.0:
            #    break


            for lapi in range(lap+conf.interval, lap+merge_step*conf.interval, conf.interval):
                print("merge lap -------", lapi)
                info['lap'] = lapi
                if not(os.path.isfile(info['particle_file'] + "_" + str(info['lap']) + ".h5")):
                    continue
                binsi, histi = create_spectra('gamma', info, conf, read_file=False, no_plots=True)
                hist += histi

            #if lap2time(lap,conf) > 2.3:
            #    break

            col = cmap(norm( lap2time(lap,conf) ))
            fit_gamma_cool(bins, hist, color=col)

            if lap == 2000:
                top_id, top_proc = top_particles_id(10, info, topgams)
            
                #Write to file 10 fastest particles
                f = open('{}/10_prtcls.txt'.format(conf.outdir), "w+")
                for i in range(10):
                    f.write("{},{},1.00\n".format(top_id[i], top_proc[i]))
                f.close()
                print("10 most energetic particles written to file")
                
    #--------------------------------------------------
    if args.lap == None:
        axleft    = 0.18
        axbottom  = 0.16
        axright   = 0.96
        axtop     = 0.83

        pos1 = axs[0].get_position()
        axwidth  = axright - axleft
        axheight = (axtop - axbottom)*0.02
        axpad = 0.01
        cax = fig.add_axes([axleft, axtop + axpad, axwidth, axheight])

        cb1 = colorbar.ColorbarBase(
                cax,
                cmap=cmap,
                norm=norm,
                orientation='horizontal',
                ticklocation='top')
        cb1.set_label(r'$t$ $(l_0/c)$')

        fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

        fname = fdir+'prtcl_spec.pdf'
        plt.savefig(fname)

    else:
        axleft    = 0.18
        axbottom  = 0.16
        axright   = 0.96
        axtop     = 0.92

        fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)
        slap = str(info['lap']).rjust(4, '0')

        fname = fdir+'prtcl_spec_{}.pdf'.format(slap)
        plt.savefig(fname)

