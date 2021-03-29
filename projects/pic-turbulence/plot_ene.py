import numpy as np
import h5py as h5
import sys, os

import matplotlib.pyplot as plt
 
from problem import Configuration_Turbulence as Configuration
from units import Units

import pytools


def read_prtcls(fname):
    f5F = h5.File(fname,'r')

    nx = f5F['Nx'][()]
    ny = f5F['Ny'][()]
    nz = f5F['Nz'][()]

    #flatten time away
    xloc = (f5F['x'][:]   ).flatten() 
    yloc = (f5F['y'][:]   ).flatten()  
    zloc = (f5F['z'][:]   ).flatten()  
    ux   = (f5F['vx'][:]  ).flatten()  
    uy   = (f5F['vy'][:]  ).flatten()  
    uz   = (f5F['vz'][:]  ).flatten()  
    wgt  = (f5F['wgt'][:] ).flatten()  
    ids  = (f5F['id'][:]  ).flatten()  
    procs= (f5F['proc'][:]).flatten() 

    #print("non zero prtcls:", np.count_nonzero(xloc))
    return xloc, yloc, zloc, ux, uy, uz, wgt, ids, procs 

# label lines packets
from math import atan2, degrees
import warnings

from matplotlib.dates import date2num, DateConverter, num2date
from matplotlib.container import ErrorbarContainer
from datetime import datetime


# Label line with line2D label data
def labelLine(line, x, label=None, align=True, offset=None, rotation=None, **kwargs):
    '''Label a single matplotlib line at position x
    Parameters
    ----------
    line : matplotlib.lines.Line
       The line holding the label
    x : number
       The location in data unit of the label
    label : string, optional
       The label to set. This is inferred from the line by default
    kwargs : dict, optional
       Optional arguments passed to ax.text
    '''
    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    mask = np.isfinite(ydata)
    if mask.sum() == 0:
        raise Exception('The line %s only contains nan!' % line)

    # Find first segment of xdata containing x
    for i, (xa, xb) in enumerate(zip(xdata[:-1], xdata[1:])):
        if min(xa, xb) <= x <= max(xa, xb):
            break
    else:
        i = len(xdata)-2
        xa = xdata[-2]
        xb = xdata[-1]
        #raise Exception('x label location is outside data range!')

    def x_to_float(x):
        """Make sure datetime values are properly converted to floats."""
        return date2num(x) if isinstance(x, datetime) else x

    xfa = x_to_float(xa)
    xfb = x_to_float(xb)
    ya = ydata[i]
    yb = ydata[i + 1]
    y = ya + (yb - ya) * (x_to_float(x) - xfa) / (xfb - xfa)

    if not(offset == None):
        y += offset

    if not (np.isfinite(ya) and np.isfinite(yb)):
        warnings.warn(("%s could not be annotated due to `nans` values. "
                       "Consider using another location via the `x` argument.") % line,
                      UserWarning)
        return

    if not label:
        label = line.get_label()

    if rotation == None:
        if align:
            # Compute the slope and label rotation
            screen_dx, screen_dy = ax.transData.transform((xfa, ya)) - ax.transData.transform((xfb, yb))
            rotation = (degrees(atan2(screen_dy, screen_dx)) + 90) % 180 - 90
            rotation=50.0
        else:
            rotation = 0


    # Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    bbox = dict(boxstyle='square', 
                pad=0.05, 
                fc=ax.get_facecolor(), 
                ec='none')

    #rotation=40.0
    ax.text(x, y, label, rotation=rotation, bbox=bbox, **kwargs)
    #ax.text(x, y, label, rotation=rotation, **kwargs)


def labelLines(lines, align=True, xvals=None, **kwargs):
    '''Label all lines with their respective legends.
    Parameters
    ----------
    lines : list of matplotlib lines
       The lines to label
    align : boolean, optional
       If True, the label will be aligned with the slope of the line
       at the location of the label. If False, they will be horizontal.
    xvals : (xfirst, xlast) or array of float, optional
       The location of the labels. If a tuple, the labels will be
       evenly spaced between xfirst and xlast (in the axis units).
    kwargs : dict, optional
       Optional arguments passed to ax.text
    '''
    ax = lines[0].axes

    labLines, labels = [], []
    handles, allLabels = ax.get_legend_handles_labels()

    all_lines = []
    for h in handles:
        if isinstance(h, ErrorbarContainer):
            all_lines.append(h.lines[0])
        else:
            all_lines.append(h)

    # Take only the lines which have labels other than the default ones
    for line in lines:
        if line in all_lines:
            label = allLabels[all_lines.index(line)]
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xvals = ax.get_xlim()  # set axis limits as annotation limits, xvals now a tuple
    if type(xvals) == tuple:
        xmin, xmax = xvals
        xscale = ax.get_xscale()
        if xscale == "log":
            xvals = np.logspace(np.log10(xmin), np.log10(xmax), len(labLines)+2)[1:-1]
        else:
            xvals = np.linspace(xmin, xmax, len(labLines)+2)[1:-1]

        if isinstance(ax.xaxis.converter, DateConverter):
            # Convert float values back to datetime in case of datetime axis
            xvals = [num2date(x).replace(tzinfo=ax.xaxis.get_units())
                     for x in xvals]

    for line, x, label in zip(labLines, xvals, labels):
        labelLine(line, x, label, align, **kwargs)

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


    axs[0].set_xlabel(r'$t/t_c$ ')
    #axs[0].set_ylabel(r'$\mathcal{E}/\mathcal{E}_{B,0}^\parallel$')
    axs[0].set_ylabel(r'$\mathcal{E}/\mathcal{E}_{F}^0$')


    tmax = 10.0
    #axs[0].set_yscale('log')
    #axs[0].set_ylim((1.0e-4, 2.0))

    axs[0].set_ylim((0.0, 1.2))
    axs[0].set_xlim((-0.5, tmax))


    pic = True

    #--------------------------------------------------
    # read and plot data
    args_cli = pytools.parse_args()
    var = args_cli.var

    #confs = [ "gam3.ini", ]
    confs = [
            args_cli.conf_filename,
            ]

    for conf_filename in confs:
        conf = Configuration(conf_filename, do_print=True)

        # restructure path to point to this dir
        fdir = conf.outdir 
        print(fdir)

        units = Units(conf)

        #--------------------------------------------------
        # read data

        flds_Epar = []
        flds_Eperp = []

        flds_Bpar = []
        flds_Bperp = []

        flds_Jpar = []
        flds_Jperp = []
        flds_time = []

        flds_tot = []

        prtcl_ene = []

        for lap in range(0, conf.Nt+1, conf.interval):
            info = {}
            info['lap'] = lap
            info['fields_file'  ] = fdir + '/flds_'+str(lap)+'.h5'

            #if lap > 2700:
            #if lap > 7500:
            #    break

            if not(os.path.isfile(info['fields_file'])):
                continue
            print(info['fields_file'])

            #print("lap {}".format(lap))
            time = units.lap2time(lap) #, do_print=True)

            #shift time with half packet step to accommodate initial positions
            time -= 0.5


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

            ene_bperp = 0.5*np.sum(bperp2) * units.norm_flds
            ene_bpar = 0.5*np.sum(bpar2) * units.norm_flds

            ene_eperp = 0.5*np.sum(eperp2) * units.norm_flds
            ene_epar = 0.5*np.sum(epar2) * units.norm_flds

            ene_b  = np.sum(bpar2 + bperp2) * units.norm_flds
            ene_e  = np.sum(eperp2 + epar2) * units.norm_flds

            flds_Bpar.append(ene_bpar)
            flds_Bperp.append(ene_bperp)

            flds_Epar.append(ene_epar)
            flds_Eperp.append(ene_eperp)

            #total energy
            tot = ene_bperp + ene_bpar + ene_epar + ene_eperp
            flds_tot.append(tot)

            #time
            flds_time.append(time)

            #end of this time step; close f5
            f5_fld.close()

            #for PIC process also particles
            #-------------------------------------------------- 
            if pic:
                # prtcl file name; check if exists
                info['prtcl_file'  ] = fdir + '/test-prtcls_'+str(lap)+'.h5'
                if not(os.path.isfile(info['fields_file'])):
                    continue
                print(info['prtcl_file'])

                xloc, yloc, zloc, ux, uy, uz, wgt, ids, procs = read_prtcls(info['prtcl_file'])
                npp = np.shape(xloc)
                #print(np.shape(xloc))

                #four vel
                u = np.sqrt(ux**2 + uy**2 + uz**2)
                gam = np.sqrt(1.0 + u**2)
                beta = u/gam 

                me = units.me
                c  = units.c
                
                #p_ene  = np.sum( gam*me*c*c * units.norm_prtcls )
                p_ene  = np.sum( gam*units.norm_prtcls )
                prtcl_ene.append(p_ene)


        flds_Bpar = np.array( flds_Bpar )
        flds_Bperp= np.array( flds_Bperp )

        flds_Epar = np.array( flds_Epar )
        flds_Eperp= np.array( flds_Eperp )

        flds_tot= np.array( flds_tot )

        prtcl_ene = np.array( prtcl_ene )

        # test magnetization
        if pic:
            print('sigmas')
            sigmas = 2.0*flds_Bpar/1.0*prtcl_ene 
            print(sigmas)

        #$ change into units of initia free energy in the packets
        ene0 = flds_Bperp[0] + flds_Eperp[0]
        flds_Bpar  /= ene0 
        flds_Bperp /= ene0 
        flds_Epar  /= ene0 
        flds_Eperp /= ene0 
        flds_tot   /= ene0 
        prtcl_ene  /= ene0 

        # show delta B par
        flds_Bpar -= flds_Bpar[0] 
        
        # delta enthalphy
        if pic:
            prtcl_ene[:] -= prtcl_ene[1]

            #prtcl_ene[:] /= prtcl_ene[1]
            #prtcl_ene[:] = np.abs(prtcl_ene - 1.0)

        # error accumulation
        #if not(conf.two_wave):
        #    flds_Bperp = np.abs(flds_Bperp / flds_Bperp[0] - 1.0)
        #    flds_Eperp = np.abs(flds_Eperp / flds_Eperp[0] - 1.0)

        flds_tot -= flds_tot[0] #difference
        flds_tot = -flds_tot #flip sign

        #flds_tot = np.abs(flds_tot)

        if False:
            print("bperp")
            print(flds_Bperp)
            print("bpar")
            print(flds_Bpar)

            print("eperp")
            print(flds_Eperp)
            print("epar")
            print(flds_Epar)
            print('etot')
            print(flds_tot)

            print('prtcl')
            print(prtcl_ene)

        #--------------------------------------------------    
        #plot field energy
        col = 'r'

        ls_E = 'solid'
        ls_B = 'solid'
        ls_Bpar = 'dashdot'

        line_fld_epar,  = axs[0].plot(flds_time, flds_Epar,  color='darkblue', linestyle='dotted', linewidth=0.7)
        line_fld_eperp, = axs[0].plot(flds_time, flds_Eperp, color='darkblue', linestyle='dashed', linewidth=0.7)

        line_fld_bpar,  = axs[0].plot(flds_time, flds_Bpar,  color='r', linestyle='dotted', linewidth=0.7)
        line_fld_bperp, = axs[0].plot(flds_time, flds_Bperp, color='r', linestyle='dashed', linewidth=0.7)

        line_tot, = axs[0].plot(flds_time, flds_tot, color='k', linestyle='solid', linewidth=0.7)

        if pic:
            line_prtcls, = axs[0].plot(flds_time, prtcl_ene, color='darkorange', linestyle='solid', linewidth=0.7)

    #simplified labels
    if True:
        fontsize = 8
        align = False
        labelLine(line_fld_bpar,  8.0, label=r'$\mathcal{E}_B^\parallel-\mathcal{E}_{B}^{\parallel,0}$',ha='center', va='center', align=align, fontsize=fontsize)
        labelLine(line_fld_epar,  8.5, label=r'$\mathcal{E}_E^\parallel$',ha='center', va='center', align=align, fontsize=fontsize)

        labelLine(line_fld_bperp, 8.0, label=r'$\mathcal{E}_B^\perp$',ha='center', va='center', align=align, fontsize=fontsize)
        labelLine(line_fld_eperp, 8.5, label=r'$\mathcal{E}_E^\perp$',ha='center', va='center', align=align, fontsize=fontsize)

        labelLine(line_tot,      10.0,   label=r'$\mathcal{E}_\nu$',ha='center', va='center', align=align, fontsize=fontsize)

        if pic:
            labelLine(line_prtcls,   10.0,   label=r'$\mathcal{E}_{p}$',ha='center', va='center', align=align, fontsize=fontsize)


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

    fname = fdir + '/ene_history.pdf'
    plt.savefig(fname)




