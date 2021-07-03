import numpy as np
import h5py as h5
import sys, os

import matplotlib.pyplot as plt
 
from configSetup import Configuration

sys.path.insert(0,'../')
from setup_tests import Configuration_Test
from setup_tests import lap2time

import pytools


if __name__ == "__main__":
    fig = plt.figure(1, figsize=(3.487, 2.35))
    
    plt.rc('font',  family='sans-serif')
    #plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.0)
    gs.update(wspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    #axs.append( plt.subplot(gs[1,0]) )
    
    for ax in axs:
        ax.minorticks_on()

    if len(axs) > 1:
        axs[0].tick_params(which='x', direction="in")
        axs[0].axes.get_xaxis().set_ticks([])

    axs[0].set_yscale('log')

    axs[0].set_xlabel(r'$t \omega_{p}$')
    axs[0].set_ylabel(r'$\frac{B^2}{8\pi n_0 \gamma_b m_e c^2}$')

    axs[0].set_ylim((1.0e-5, 1.0))
    axs[0].set_xlim((0, 100.0))


    #--------------------------------------------------
    # read and plot data

    confs = [
            "gam3.ini",
            "gam10.ini",
            "gam30.ini",
            "gam100.ini",
            ]


    for ii, conf_filename in enumerate(confs):
        conf = Configuration_Test(conf_filename)

        # restructure path to point to this dir
        fdir = conf.outdir 
        fdir = fdir.split('/')[1] + '/'
        print(fdir)

        col = "C"+str(ii)

        #--------------------------------------------------
        # theoretical growth rate
        gam = conf.gamma
        beta = np.sqrt(1.0 - 1.0/gam**2)

        delta = 1.0*beta*np.sqrt(2/gam)
        print("growth rate:", delta)

        tt = np.linspace(5.0, 100.0, 100)
        yy = np.exp(delta*tt)*2.0e-6

        print(tt)
        print(yy)
        axs[0].plot(tt, yy, linestyle="--", color='k', alpha=0.9)
        if ii == 0:
            axs[0].text(10.0, 0.1, r'$\gamma_{\mathrm{b}} = 3$', rotation=80, fontsize = 7, ha='center', va='center')
        if ii == 1:
            axs[0].text(23.0, 0.3, r'$\gamma_{\mathrm{b}} = 10$', rotation=72, fontsize = 7, ha='center', va='center')
        if ii == 2:
            axs[0].text(42.0, 0.3, r'$\gamma_{\mathrm{b}} = 30$', rotation=60, fontsize = 7, ha='center', va='center')
        if ii == 3:
            axs[0].text(50.0, 0.004, r'$\gamma_{\mathrm{b}} = 100$', rotation=45, fontsize = 7, ha='center', va='center')


        #--------------------------------------------------
        # normalization and units
        qe = np.abs(conf.qe)
        me = np.abs(conf.me)*qe
        c  = conf.cfl 
        print("me: {} qe: {} c: {}".format(me, qe, c))

        Lx, Ly, Lz = conf.NxMesh*conf.Nx, conf.NyMesh*conf.Ny, conf.NzMesh*conf.Nz
        norm_prtcl = Lx*Ly*Lz*conf.ppc*2.0 #pair-plasma so multiply by 2
        norm_prtcl *= conf.gamma*me*c**2 #kinetic energy in flow

        norm_flds = 1.0
        norm_flds /= 2.0*(1.0/conf.stride)*(1.0/conf.stride) #divide by volume of simulation domain

        #3D switch
        if Lz > 1:
            norm_flds /= (1.0/conf.stride)

        # fields into units of particle enthalphy density
        norm_flds /= norm_prtcl


        #--------------------------------------------------
        # read data

        flds_E = []
        flds_B = []
        flds_Bpar = []
        flds_time = []

        for lap in range(0, conf.Nt+1, conf.interval):
            info = {}
            info['lap'] = lap
            info['fields_file'  ] = fdir + 'flds_'+str(lap)+'.h5'

            print(info['fields_file'])
            if not(os.path.isfile(info['fields_file'])):
                continue
            
            if lap > 800:
                continue

            #print("lap {}".format(lap))
            time = lap2time(lap, conf) #, do_print=True)

            # read file
            f5_fld = h5.File(info['fields_file'],'r')
            bx  = pytools.read_h5_array(f5_fld, 'bx')
            by  = pytools.read_h5_array(f5_fld, 'by')
            bz  = pytools.read_h5_array(f5_fld, 'bz')

            jx  = pytools.read_h5_array(f5_fld, "jx")
            jy  = pytools.read_h5_array(f5_fld, "jy")
            jz  = pytools.read_h5_array(f5_fld, "jz")

            ex  = pytools.read_h5_array(f5_fld, "ex")
            ey  = pytools.read_h5_array(f5_fld, "ey")
            ez  = pytools.read_h5_array(f5_fld, "ez")

            bperp2 = (bx*bx + by*by) #B_perp^2
            bpar2  = (bz*bz) #B_par^2

            eperp2 = (ex*ex + ey*ey) #E_perp^2
            epar2  = (ez*ez) #E_par^2

            ene_bperp = 0.5*np.sum(bperp2) * norm_flds
            ene_bpar = 0.5*np.sum(bpar2) * norm_flds

            ene_b  = np.sum(bpar2 + bperp2) * norm_flds
            ene_e  = np.sum(eperp2 + epar2) * norm_flds


            flds_B.append(ene_b)
            flds_E.append(ene_e)

            flds_time.append(time)

            #end of this time step; close f5
            f5_fld.close()

        #--------------------------------------------------    
        #plot field energy
        #col = 'r'

        ls_E = 'solid'
        ls_B = 'solid'
        ls_Bpar = 'dashdot'


        print(flds_time)
        print(flds_B)

        #line_fld_e, = axs[0].plot(flds_time, flds_E,  color='darkblue', linestyle=ls_E, linewidth=0.7)
        line_fld_b, = axs[0].plot(flds_time, flds_B,  color=col, linestyle=ls_B, linewidth=1.0)


    #--------------------------------------------------    
    # close and save fig
    axleft    = 0.18
    axbottom  = 0.16
    axright   = 0.96
    axtop     = 0.93

    pos1 = axs[0].get_position()
    axwidth  = axright - axleft
    axheight = (axtop - axbottom)*0.02
    axpad = 0.01

    #cax = fig.add_axes([axleft, axtop + axpad, axwidth, axheight])
    #cb1 = colorbar.ColorbarBase(
    #        cax,
    #        cmap=cmap,
    #        norm=norm,
    #        orientation='horizontal',
    #        ticklocation='top')
    #cb1.set_label(r'$t$ $(l_0/c)$')

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    fname = 'weibel_v2.pdf'
    plt.savefig(fname)




