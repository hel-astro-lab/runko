import numpy as np
import h5py as h5
import sys, os
import argparse

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import colorbar
from scipy.optimize import minimize


from parser import parse_input
from configSetup import Configuration
from problem import Configuration_Problem
from testprtcls import TestParticles



def lap2time(lap, conf, do_print=False):

    # from relativistic pair plasma oscillation freq to single species classical
    omp = conf.omp*np.sqrt(conf.gamma/2)

    return lap*omp


#from labellines import labelLine, labelLines
#from plot_spec_fit import draw_label

def read_var(f5F, var, conf):
    try:
        val = f5F[var][:,:,0]
    except:
        nx = f5F['Nx'][()]
        ny = f5F['Ny'][()]
        nz = f5F['Nz'][()]
        #print("reshaping 1D array into multiD with {} {} {}".format(nx,ny,nz))

        #print(nx,ny,nz)
        val = f5F[var][:]
        val = np.reshape(val, (nx, ny, nz))
        val = val[:,:,0]

    return val.astype(float)


if __name__ == "__main__":
    fig = plt.figure(1, figsize=(3.487, 2.5))
    #fig = plt.figure(1, figsize=(3.487, 4.0))
    
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

    #--------------------------------------------------
    if False:
        left, bottom, width, height = [0.65, 0.67, 0.3, 0.15]
        ax2 = fig.add_axes([left, bottom, width, height])
        #ax2.set_yscale('log')
        #ax2.set_xscale('log')
        ax2.set_xlabel(r"$t \omega_p$", fontsize=5)
        ax2.set_ylabel(r"$\theta$", fontsize=5)
        ax2.xaxis.set_tick_params(labelsize=4)
        ax2.yaxis.set_tick_params(labelsize=4)
        ax2.tick_params(which='both', direction="in")


    #--------------------------------------------------
    # command line driven version
    conf, fdir, args = parse_input()
    conf = Configuration_Problem(args.conf_filename)

    fdir += '/'
    print("processing {}".format(fdir))


    #--------------------------------------------------
    # set up plot
    tmin = 0.0
    tmax = 20.0

    ax.set_yscale('log')
    #ax.set_xscale('log')

    #ax.set_xlim((tmin, tmax))
    #ax.set_ylim((0.0, 1.0))
    
    ax.set_xlabel(r'$t \omega_p$')
    #ax.set_ylabel(r"$E / E_{\mathrm{mag},\perp}(t=0)$")


    # theoretical growth rate
    gam = conf.gamma
    beta = np.sqrt(1.0 - 1.0/gam**2)

    delta = 1.0*beta*np.sqrt(2/gam)
    #delta = np.sqrt(3.0/gam)/2
    print("Growth rate:", delta)

    tt = np.linspace(0.0, 70., 100)
    yy = np.exp(delta*tt) * 5.0e-6
    ax.plot(tt, yy, "k--")


    #--------------------------------------------------
    # normalization and units
    qe = np.abs(conf.qe)
    me = np.abs(conf.me)*qe
    c  = conf.cfl 
    print("me: {} qe: {} c: {}".format(me, qe, c))

    #--------------------------------------------------
    #grid points
    Lx, Ly, Lz = conf.NxMesh*conf.Nx, conf.NyMesh*conf.Ny, conf.NzMesh*conf.Nz
    #Benerg = 0.5*Lx*Ly*Lz*conf.binit**2
    #Benerg = 1.0

    #--------------------------------------------------
    #particle number
    #n_test_prtcls = 999424.0
    #n_test_prtcls = 102400

    norm_prtcl = Lx*Ly*Lz*conf.ppc*2.0 #pair-plasma so multiply by 2
    print("total number of particles: ", norm_prtcl)
    norm_prtcl *= conf.gamma*me*c**2 #kinetic energy in flow

    #norm_prtcl /= n_test_prtcls #1.0e6 #fixed number of test particles used
    print("prtcl/test_prtcl ratio: ", norm_prtcl)
    #norm_prtcl /= Benerg


    #--------------------------------------------------
    #field normalization; NEW with every stride:th value stored 
    #norm_flds = 2.0/conf.binit**2 #mean initial energy density
    #norm_flds = 2.0
    #norm_flds /= (Lx/conf.stride)*(Ly/conf.stride) #divide by volume of simulation domain
    norm_flds = 1.0
    norm_flds /= 2.0*(1.0/conf.stride)*(1.0/conf.stride) #divide by volume of simulation domain

    #3D switch
    if Lz > 1:
        norm_flds /= (1.0/conf.stride)
    norm_flds /= norm_prtcl


    # for actual particles use sum of particles multiplied by gamma
    n_test_prtcls = conf.n_test_prtcls
    norm_prtcl = 1.0/(conf.gamma*n_test_prtcls*me*c**2)


    #--------------------------------------------------
    print("====normalizations:=====")
    print("norm_prtcl: ", norm_prtcl)
    print("norm_flds: ",  norm_flds)

    print("-------checking units--------")
    #units = Units(conf)
    #print(units.norm_flds/norm_flds)
    #print(units.norm_prtcls/(norm_prtcl*me*c*c))

    break_lap = 400

    #--------------------------------------------------
    # visualization 
    #norm = Normalize(vmin=0.0, vmax=10.)
    #cmap = get_cmap('Spectral')


    prtcl_th   = []
    prtcl_nth  = []
    prtcl_time = []

    flds_E = []
    flds_B = []
    flds_Bpar = []
    flds_time = []

    #--------------------------------------------------
    # prtcl energetics

    #--------------------------------------------------
    # save to/load from file
    prtcls = TestParticles()
    if False:
        ####################
        # read prtcl data
        fname_P = "test-prtcls"

        i = 0
        files_P = []
        for lap in range(0, conf.Nt+1, conf.interval):
            #if lap > 5:
            #    break

            files_P.append( fdir+fname_P)

            info = {}
            info['lap'] = lap
            info['particle_file'] = fdir+fname_P
            prtcls.add_file(info)
            i += 1

        #print("n prtcls", len(prtcls.prtcls)) 
        #print(prtcls.times)
        prtcls.times = np.array(prtcls.times)

        ####################
        # build gamma array

        gams = np.zeros((len(prtcls.prtcls), len(prtcls.times)))

        #print("gams size:", np.shape(gams))

        for ip,key in enumerate(prtcls.prtcls):
            tslices = prtcls.prtcls[key]

            #print(tslices)
            #g = tslices[7] #get gamma as a function of time

            g = np.zeros(len(tslices))
            for j,prtcl in enumerate(tslices):
                #ux = prtcl[3]
                #uy = prtcl[4]
                #uz = prtcl[5]
                #gam = np.sqrt(1.0 + ux*ux + uy*uy + uz*uz)
                g[j] = prtcl[7]
            gams[ip,:] = g

        #--------------------------------------------------
        # save to/load from file

        #write to file
        #if True:
        f5 = h5.File(fdir + 'prtcl_data.h5', 'w')
        dset1  = f5.create_dataset("gams",  data=gams)
        dset3  = f5.create_dataset("times", data=prtcls.times)
        f5.close()

    #read from file
    else:
        f5 = h5.File(fdir + 'prtcl_data.h5', 'r')
        gams = f5['gams'][:,:]
        prtcls.times = f5['times'][:]
        f5.close()


    ####################
    # build energy arrays from individual prtcl data
    if True:
        print("shape of th", np.shape(gams))

        for it, lap in enumerate(prtcls.times):
            if lap > break_lap:
                break

            gam1 = gams[:,it]
            #gam1 = gam1[np.nonzero(gam1)]
            beta1 = np.sqrt(1.0 - 1.0/gam1**2)

            if lap % 10 == 0:
                #print("mean gamma:", np.mean(gam1))
                t = lap2time(lap, conf, do_print=True)
                #print("n(th): {} n(nth): {} n_nth/n_tot: {}".format(
                #len(gam1), len(gam2), len(gam2)/(len(gam1)+len(gam2))))
            else:
                t = lap2time(lap, conf)

            ###################################################################
            #particle energy
            ene_th  = np.sum( gam1*me*c*c * norm_prtcl )
            prtcl_th.append(ene_th)
            prtcl_time.append(t)



    #--------------------------------------------------
    # field energetics

    #read field files and calculate energy stored in magnetic and electric fields
    if True:

        it = -1
        for lap in range(0, conf.Nt+1, conf.interval):
            if lap > break_lap:
                break
            it += 1

            info = {}
            info['lap'] = lap
            info['fields_file'  ] = fdir + 'flds_'+str(lap)+'.h5'

            if not(os.path.isfile(info['fields_file'])):
                continue

            #print("lap {}".format(lap))
            if lap % 20 == 0:
                time = lap2time(lap, conf, do_print=True)
            else:
                time = lap2time(lap, conf)

            # read file
            f5_fld = h5.File(info['fields_file'],'r')
            bx  = read_var(f5_fld, 'bx', conf)
            by  = read_var(f5_fld, 'by', conf)
            bz  = read_var(f5_fld, 'bz', conf)

            jx  = read_var(f5_fld, "jx", conf)
            jy  = read_var(f5_fld, "jy", conf)
            jz  = read_var(f5_fld, "jz", conf)

            ex  = read_var(f5_fld, "ex", conf)
            ey  = read_var(f5_fld, "ey", conf)
            ez  = read_var(f5_fld, "ez", conf)

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
    # plotting options

    col = 'k'
    ls_th = 'solid'
    ls_nth = 'dashed'

    ls_E = 'solid'
    ls_B = 'solid'
    ls_Bpar = 'dashdot'

    #col = cmap(norm( lap2time(it*conf.interval,conf) ))

    #--------------------------------------------------    
    #particle kinetic + rest mass
    line_prtcl_th,  = axs[0].plot(prtcl_time, prtcl_th,   color='k', linestyle=ls_th,  linewidth=0.7)


    #--------------------------------------------------    
    #field energy
    col = 'r'
    line_fld_e, = axs[0].plot(flds_time, flds_E,  color='darkblue', linestyle=ls_E, linewidth=0.7)
    line_fld_b, = axs[0].plot(flds_time, flds_B,  color='red',      linestyle=ls_B, linewidth=0.7)
    #line_fld_bp,= axs[0].plot(flds_time, flds_Bpar,  color='red',   linestyle=ls_Bpar, linewidth=0.7)



    #simplified labels
    if False:
        fontsize = 6
        align = True
        labelLine(line_prtcl_th,  9.0, label=r'prtcl (th)',        ha='center', va='center', align=align, fontsize=fontsize)
        labelLine(line_prtcl_nth, 7.0, label=r'prtcl (nth)',       ha='center', va='center', align=align, fontsize=fontsize)

        labelLine(line_rad_th,  11.0,  label=r'rad (th)',          ha='center', va='center', align=align, fontsize=fontsize)
        labelLine(line_rad_nth,  9.0,  label=r'rad (nth)',         ha='center', va='center', align=align, fontsize=fontsize)

        labelLine(line_fld_b,    4.0,  label=r'mag ($\perp$)',     ha='center', va='center', align=align, fontsize=fontsize)
        labelLine(line_fld_e,    2.5,  label=r'elec',              ha='center', va='center', align=align, fontsize=fontsize)
        #labelLine(line_fld_bp,  15.0,  label=r'mag ($\parallel$)', ha='center', va='center', align=align, fontsize=fontsize)



    #--------------------------------------------------    
    #sum of everything
    ene_0 = flds_B[0] + flds_E[0] + prtcl_th[0]  #total energy content at the beginning


    #real sigma in simulation
    #print()
    #print("real sigma at t=0:", 2.0*flds_B[0] / (prtcl_th[0] + prtcl_nth[0]))
    #print()

    #axs[0].axhline(ene_0, color='k', linestyle='dotted', linewidth=1.0, alpha=0.6)
    #axs[0].axhline(1.0+1.0/conf.sigma, color='b', linestyle='dotted', linewidth=1.0, alpha=0.6)

    
    flds_time = np.array(flds_time)
    prtcl_time = np.array(prtcl_time)

    ene_sum = []
    sum_time = []

    lost_ene = []

    for ip, t in enumerate(prtcl_time):
        ifld = np.where(flds_time == t)[0]
        if len(ifld) == 0:
            continue
        ifld = ifld[0]

        esum  = prtcl_th[ip] 
        esum += flds_B[ifld] + flds_E[ifld] 

        #lost energy
        loste = ene_0 - esum 
        lost_ene.append(loste)

        #esum += rad_th[ip] + rad_nth[ip] 
        ene_sum.append(esum)
        sum_time.append(t)

        #print("t: {} \t lost E: {}".format(t, loste))


    line_tot, = axs[0].plot(sum_time, ene_sum, color='k', linestyle='dashed', linewidth=1.3, alpha=0.9)


    #--------------------------------------------------    
    # write data points to file
    if True:
        f5 = h5.File(fdir + 'ene.h5', 'w')

        #prtcl_time
        f5.create_dataset("prtcl_th",  data=prtcl_th)
        f5.create_dataset("prtcl_time",data=prtcl_time)

        f5.create_dataset("flds_E",    data=flds_E)
        f5.create_dataset("flds_B",    data=flds_B)
        f5.create_dataset("flds_time", data=flds_time)

        f5.close()


    #--------------------------------------------------    
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

    fname = fdir+'fld_evolution.pdf'
    plt.savefig(fname)



