import numpy as np
import matplotlib.pyplot as plt

import h5py as h5

import sys, os
#from parser import parse_input
import argparse
import matplotlib
import pytools  # runko python tools

from init_problem import Configuration_Turbulence as Configuration



#from utils_prtcl_spec import default_values
#from utils_prtcl_spec import default_turbulence_values
#from utils_prtcl_spec import create_spectra


from utils_analysis import lap2time
from utils_analysis import create_pmom_spectra



if __name__ == "__main__":

    fig = plt.figure(1, figsize=(3.487, 5.5))

    plt.rc('font',  family='sans-serif')
    #plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(2, 1)
    gs.update(hspace = 0.4)
    #gs.update(wspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    axs.append( plt.subplot(gs[1,0]) )
    #axs.append( plt.subplot(gs[2,0]) )
    
    for ax in axs:
        ax.minorticks_on()
        ax.set_xlim((1e-1, 1e3))

        ax.set_yscale('log')
        ax.set_xscale('log')

    #axs[0].set_xlim((1.0, 1000.0))
    #axs[1].set_xlim((1.0, 1000.0))
    #axs[2].set_xlim((1.0, 1000.0))
                           
    axs[0].set_xlabel(r"$p$ ($m_e c$)")
    axs[0].set_ylabel(r"$p \mathrm{d} n_-/\mathrm{d} p$")

    axs[1].set_xlabel(r"$p$ ($m_e c$)")
    axs[1].set_ylabel(r"$p \mathrm{d} n_+/\mathrm{d} p$")

    #axs[2].set_xlabel(r"$x$ ($m_e c^2$)")
    #axs[2].set_ylabel(r"$x \mathrm{d} n_x/\mathrm{d} x$")

    #--------------------------------------------------
    # command line driven version
    #conf, fdir, args = parse_input()
    #fdir += '/'

    # parse command line arguments
    do_print = True
    args = pytools.parse_args()
    conf = Configuration(args.conf_filename, do_print=do_print)
    fdir = conf.outdir + "/"


    print("processing {}".format(fdir))
    fname_P = "test-prtcls"

    sp_args = {}
    sp_args['xmin'] = 1e-1
    sp_args['xmax'] = 1e+3
    sp_args['nbins'] = 64


    # if lap is defined, only process that one individual round
    if not(args.lap == None):
        laptime = lap2time(args.lap, conf)

        for ispc in [0,1]:
            if ispc == 0:
                fname = fdir + fname_P + "_" + str(args.lap) + ".h5"
            else:
                fname = fdir + fname_P + "-" + str(ispc) + "_" + str(args.lap) + ".h5"
            if not(os.path.isfile(fname)): continue

            print('fname', fname, 't:', laptime, 'isp:', ispc)
            bins, hist, histx, histy, histz = create_pmom_spectra(fname, sp_args)

            #p2_dn_dp = bins*hist # p^2 dn/dp = p dn/dlog(p)
            p_dn_dp = hist # p dn/dp = dn/dlog(p)
                
            axs[ispc].plot(bins, p_dn_dp, color='k', linewidth=0.8, linestyle='solid', alpha=0.9)


    # else process every file that there is
    else:
        tmax = 5.0
        norm = matplotlib.colors.Normalize(vmin=0.0, vmax=tmax)
        #cmap = matplotlib.colormaps['Spectral']
        cmap = matplotlib.colormaps['turbo_r']

        files_P = []

        merge_step = 1
        for lap in range(0, conf.Nt+1, merge_step*conf.interval):
            laptime = lap2time(lap, conf)

            col = cmap(norm(lap2time(lap, conf)))

            if laptime > tmax:
                continue

            for ispc in [0,1]:

                if ispc == 0:
                    fname = fdir + fname_P + "_" + str(lap) + ".h5"
                else:
                    fname = fdir + fname_P + "-" + str(ispc) + "_" + str(lap) + ".h5"

                if not(os.path.isfile(fname)): continue

                print('fname', fname, 't:', laptime, 'isp:', ispc)
                bins, hist, histx, histy, histz = create_pmom_spectra(fname, sp_args)

                #p2_dn_dp = bins*hist # p^2 dn/dp = p dn/dlog(p)
                p_dn_dp = hist # p^2 dn/dp = p dn/dlog(p)
                    
                axs[ispc].plot(bins, p_dn_dp, color=col, linewidth=0.8, linestyle='solid', alpha=0.9)
                        


    #--------------------------------------------------
    if args.lap == None:
        axleft    = 0.18
        axbottom  = 0.10
        axright   = 0.96
        axtop     = 0.92

        pos1 = axs[0].get_position()
        axwidth  = axright - axleft
        axheight = (axtop - axbottom)*0.02
        axpad = 0.01
        cax = fig.add_axes([axleft, axtop + axpad, axwidth, axheight])

        cb1 = matplotlib.colorbar.ColorbarBase(
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
        axbottom  = 0.10
        axright   = 0.96
        axtop     = 0.92

        fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)
        slap = str(args.lap).rjust(4, '0')

        fname = fdir+'prtcl_spec_{}.pdf'.format(slap)
        plt.savefig(fname)

