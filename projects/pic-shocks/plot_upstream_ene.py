import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import sys, os 

from init_problem import Configuration_Shocks as Configuration_Problem

import pytools
from pytools import read_h5_array
from pytools.cli import parse_args

from matplotlib.colors import Normalize
#from matplotlib.cm import get_cmap
from matplotlib import colorbar

from labellines import labelLine, labelLines

from init_problem import lap2time


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
        ax2.set_xlabel(r"$t$ $(l_0/c)$", fontsize=5)
        ax2.set_ylabel(r"$\theta$", fontsize=5)
        ax2.xaxis.set_tick_params(labelsize=4)
        ax2.yaxis.set_tick_params(labelsize=4)
        ax2.tick_params(which='both', direction="in")


    #--------------------------------------------------
    # command line driven version
    args = pytools.parse_args()
    conf = Configuration_Problem(args.conf_filename)
    fdir = conf.outdir + '/'

    print("processing {}".format(fdir))



    #--------------------------------------------------
    # units
    comp = conf.c_omp
    b0 = conf.binit
    beta = conf.beta
    
    print('comp', comp)
    print('b0', b0)
    print('beta', beta)

    #--------------------------------------------------
    # set up plot
    tmin = 0.0
    tmax = 5e3

    ax.set_yscale('log')
    #ax.set_xscale('log')

    ax.set_xlim((tmin, tmax))
    #ax.set_ylim((1e-9, 1e1))
    ax.set_ylim((1e-5, 1e0))
                           
    ax.set_xlabel(r'$t \omega_p$')
    ax.set_ylabel(r"$dB_i^2/B_0^2$")


    #--------------------------------------------------
    break_lap = 300000


    #--------------------------------------------------
    # visualization 
    norm = Normalize(vmin=0.0, vmax=tmax)
    #cmap = get_cmap('Spectral')

    
    tt  = []

    xib = []
    xib2 = []

    bxd = []
    byd = []
    bzd = []

    exd = []
    eyd = []
    ezd = []

    poynt = []

    prev_Sx = 1.0e-2
    for lap in range(0, conf.Nt+1, conf.interval):
        time = lap2time(lap,conf)
        #col = cmap(norm( time ))

        fname = fdir + "shock_" + str(lap) + ".h5"
        if not(os.path.isfile(fname)):
            continue
        print(time, fname)

        # read data
        f5F = h5.File(fname,'r')

        #limit to a box ahead of shock
        #xshift1 = int((-conf.box_shift + 10.0*comp)/conf.box_stride)
        #xshift2 = int((-conf.box_shift + 10.0*comp + 40.0*comp)/conf.box_stride)

        # Plotinkov & Sironi range
        xshift1 = int((-conf.box_shift + 5.0*comp)/conf.box_stride)
        xshift2 = int((-conf.box_shift + 5.0*comp + 25.0*comp)/conf.box_stride)
        #print(xshift)

        rho=read_h5_array(f5F, 'rho')[xshift1:xshift2, :, :]
        ex = read_h5_array(f5F, 'ex')[xshift1:xshift2, :, :] + conf.ex_ext
        ey = read_h5_array(f5F, 'ey')[xshift1:xshift2, :, :] + conf.ey_ext
        ez = read_h5_array(f5F, 'ez')[xshift1:xshift2, :, :] + conf.ez_ext

        bx = read_h5_array(f5F, 'bx')[xshift1:xshift2, :, :] + conf.bx_ext
        by = read_h5_array(f5F, 'by')[xshift1:xshift2, :, :] + conf.by_ext
        bz = read_h5_array(f5F, 'bz')[xshift1:xshift2, :, :] + conf.bz_ext

        jx = read_h5_array(f5F, 'jx')[xshift1:xshift2, :, :]
        jy = read_h5_array(f5F, 'jy')[xshift1:xshift2, :, :]
        jz = read_h5_array(f5F, 'jz')[xshift1:xshift2, :, :]

        # projected density 
        dens = f5F['dens_profile']
        sloc = int( f5F.attrs['shock_loc'][()] )

        f5F.close()

        nx,ny,nz = np.shape(rho)
        #print('box_size', nx,ny,nz)
        #print('box_size in comp' , nx/comp, ny/comp, nz/comp)

        #-------------------------------------------------- 
        # filter out some values
        #if False and conf.bperp > 0.0:
        #    if conf.bperp > 0.0:
        #        overlapind = np.where(bz/conf.binit > 1.6)
        #    elif conf.bplan > 0.0:
        #        overlapind = np.where(by/conf.binit > 1.6)

        #    ex[overlapind] /= 2.0
        #    ey[overlapind] /= 2.0
        #    ez[overlapind] /= 2.0

        #    bx[overlapind] /= 2.0
        #    by[overlapind] /= 2.0
        #    bz[overlapind] /= 2.0

        #-------------------------------------------------- 
        # analyze

        tt.append(time)

        # out-of-the-plane configuration
        #if conf.bplan == 0.0:
        if False:
            dbz = bz-b0
            dby = by
            dey = ey - beta*b0
            xib.append( np.mean( dbz*dey/b0**2 ) )

            bxd.append( np.mean( bx**2/b0**2 ))
            byd.append( np.mean( dby**2/b0**2 )) # O mode
            bzd.append( np.mean( dbz**2/b0**2 ) ) #X mode

            dez = ez
            poynt.append( np.mean(dey*dbz + dez*dby)/b0**2 )

            #bzdp.append( np.mean( (bz-b0)**2/b0**2 ) ) #1D X mode proxy
        else:
        #general configuration

            # background guide fields in all orientations
            bx0 = b0*conf.bpar
            by0 = b0*conf.bplan
            bz0 = b0*conf.bperp

            # fluctuating B fields
            dbx = bx - bx0
            dby = by - by0
            dbz = bz - bz0

            # motional electric field
            # Emot = v x B= (since v = vx)
            #              xdir ( vy*bz - vz*by )
            #              ydir (-vx*bz + vz*bx )
            #              zdir ( vx*by - vy*bx )
            #             =
            #              xdir (   0   -   0   )
            #              ydir (-vx*bz +   0   )
            #              zdir ( vx*by -   0   )
            # so vector components are:
            # Emot_x = 0.0
            # Emot_y = -vx*bz
            # Emot_z = +vx*by

            ex0 = 0.0
            ey0 = -beta*bz0
            ez0 = +beta*by0

            # fluctuating electric fields; TODO: check - sign
            dex = ex - ex0
            dey = ey - ey0
            dez = ez - ez0

            #print('minmax', np.min(dbz), np.max(dbz))
            #print('mean rho', np.mean(rho))

            xib.append(  np.mean(+dbz*dey/b0**2 ) )
            xib2.append( np.mean(-dby*dez/b0**2 ) )

            bxd.append( np.mean( (dbx**2)/b0**2 ))
            byd.append( np.mean( (dby**2)/b0**2 )) # O mode
            bzd.append( np.mean( (dbz**2)/b0**2 )) # X mode

            exd.append( np.mean( (dex**2)/b0**2 ))
            eyd.append( np.mean( (dey**2)/b0**2 )) 
            ezd.append( np.mean( (dez**2)/b0**2 )) 

            # poynting flux to x direction; TODO; check - sign
            Sx = np.mean(dey*dbz + dez*dby)/b0**2 
            if Sx > 0:
                poynt.append( Sx ) 
                prev_Sx = Sx
            else:
                poynt.append( prev_Sx ) 


    #--------------------------------------------------    
    tt    = np.array(tt)
    xib   = np.array(xib)
    xib2  = np.array(xib2)
    bxd   = np.array(bxd)
    byd   = np.array(byd)
    bzd   = np.array(bzd)
    exd   = np.array(exd)
    eyd   = np.array(eyd)
    ezd   = np.array(ezd)
    poynt = np.array(poynt)

    #--------------------------------------------------    
    # plotting options

    print('xib', xib)
    print('xib2', xib2)
    print('poynt', poynt)
    print()
    print('dbx', bxd)
    print('dby', byd)
    print('dbz', bzd)
    print()
    print('dex', exd)
    print('dey', eyd)
    print('dez', ezd)

    #--------------------------------------------------    
    #particle kinetic + rest mass
    #line_prtcl_th,  = axs[0].plot(prtcl_time, prtcl_th,   color='k', linestyle=ls_th,  linewidth=0.7)

    lw = 0.4
    alpha = 0.6
    line_xib, = axs[0].plot(tt, xib,  color='b', linestyle='dotted', alpha=alpha, linewidth=1.0, label=r"$\delta B_z \delta E_y$")
    line_xib2,= axs[0].plot(tt, xib2, color='y', linestyle='dotted', alpha=alpha, linewidth=1.0, label=r"$\delta B_y \delta E_z$")

    line_poy, = axs[0].plot(tt, poynt,color='cyan', linestyle='-',   alpha=alpha, linewidth=1.0, label=r"$S_x$")

    line_bx,  = axs[0].plot(tt, bxd,  color='r', linestyle='-',   alpha=alpha, linewidth=lw , label=r"$\delta B_x^2$")
    line_by,  = axs[0].plot(tt, byd,  color='y', linestyle='-',   alpha=alpha, linewidth=lw , label=r"$\delta B_y^2$")
    line_bz,  = axs[0].plot(tt, bzd,  color='b', linestyle='-',   alpha=alpha, linewidth=lw , label=r"$\delta B_z^2$")

    line_ex,  = axs[0].plot(tt, exd,  color='r', linestyle='dashed',   alpha=alpha, linewidth=lw , label=r"$\delta E_x^2$")
    line_ey,  = axs[0].plot(tt, eyd,  color='y', linestyle='dashed',   alpha=alpha, linewidth=lw , label=r"$\delta E_y^2$")
    line_ez,  = axs[0].plot(tt, ezd,  color='b', linestyle='dashed',   alpha=alpha, linewidth=lw , label=r"$\delta E_z^2$")

    if conf.bplan > 0:
        line_ox,  = axs[0].plot(tt, np.array(bzd)/np.array(byd),  color='green', alpha=alpha, linestyle='-', linewidth=1.0, label=r"$B_y^2/B_z^2$")
    else:
        line_ox,  = axs[0].plot(tt, np.array(byd)/np.array(bzd),  color='green', alpha=alpha, linestyle='-', linewidth=1.0, label=r"$B_y^2/B_z^2$")


    #simplified labels
    if False:
        fontsize = 6
        align = False
        labelLine(line_xib,  8000.0, label=r'$\xi_B$', ha='center', va='center', align=align, fontsize=fontsize)

        labelLine(line_bx,   9000.0, label=r'$\delta B_x/B_0$', ha='center', va='center', align=align, fontsize=fontsize)
        labelLine(line_by,   9000.0, label=r'$\delta B_y/B_0$', ha='center', va='center', align=align, fontsize=fontsize)
        labelLine(line_bz,   9000.0, label=r'$\delta B_z/B_0$', ha='center', va='center', align=align, fontsize=fontsize)

        labelLine(line_ox,   9000.0, label=r'ratio O/X', ha='center', va='center', align=align, fontsize=fontsize)


    # normal legends
    if True:
        Leg1 = plt.legend(handles=[
                line_xib,
                line_xib2,
                line_poy,
                line_bx,
                line_by,
                line_bz,
                line_ex,
                line_ey,
                line_ez,
                line_ox,
                ], 
                loc='lower right',
                #mode="expand", 
                borderaxespad=1, 
                ncol=2,
                fontsize=4,
                )
        ax = plt.gca().add_artist(Leg1)


    #--------------------------------------------------    
    # write data points to file
    if True:
        f5 = h5.File(fdir + 'upstream_ene_5_25.h5', 'w')

        f5.create_dataset("tt"   , data=tt)
        f5.create_dataset("xib"  , data=xib)
        f5.create_dataset("xib2" , data=xib2)
        f5.create_dataset("bxd"  , data=bxd)
        f5.create_dataset("byd"  , data=byd)
        f5.create_dataset("bzd"  , data=bzd)
        f5.create_dataset("exd"  , data=exd)
        f5.create_dataset("eyd"  , data=eyd)
        f5.create_dataset("ezd"  , data=ezd)
        f5.create_dataset("poynt", data=poynt)

        f5.close()

    #--------------------------------------------------    
    axleft    = 0.18
    axbottom  = 0.18
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

    fname = fdir+'ene_upstream_evolution_5_25.pdf'
    plt.savefig(fname)



