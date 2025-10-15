# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os
import matplotlib.ticker as ticker

from init_problem import Configuration_Shocks as Configuration_Problem

import pytools
from pytools import read_h5_array
from pytools.cli import parse_args

from init_problem import lap2time

# trick to make nice colorbars
# see http://joseph-long.com/writing/colorbars/
def colorbar(mappable, 
        loc="right", 
        orientation="vertical", 
        size="1%", 
        #pad=0.05, 
        pad=0.1, 
        ticklocation='right'):
        #loc="top", 
        #orientation="horizontal", 
        #size="8%", 
        #pad=0.03, 
        #ticklocation='top'):

    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(loc, size=size, pad=pad)
    return cax, fig.colorbar(mappable, cax=cax, orientation=orientation, ticklocation=ticklocation)


default_values = {
    'cmap':"viridis",
    'vmin': None,
    'vmax': None,
    'clip': None,
    'aspect':1.0, #picture aspect ratio
    'crop_ratio':8.0, # crop ratio
    'vsymmetric':None,
    'winsorize_min':0.005,
    'winsorize_max':0.005,
    'title':'',
    'log':False,
    'derived':False,
    'streamlines':False,
    'file':'fields_file',
}

default_problem_values = {
        'rho': {'title': r"$n/n_0$",
                'vmin': 0.0,
                'vmax': 4.0,
                },
        'logrho': {'title': r"$\log_{10} n/n_0$",
                'vmin': -1.0,
                'vmax':  1.0,
                'derived':True,
                'log':True,
                },
        'jz': {'title': r"$J_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
               'vmin': -0.001,
               'vmax':  0.001,
                },
        'jx': {'title': r"$J_x$",
               'cmap': "RdBu",
               'vsymmetric':True,
               #'vmin': -2.0000,
               #'vmax':  2.0000,
               'vmin': -0.5000,
               'vmax':  0.5000,
                },
        'jy': {'title': r"$J_y$",
               'cmap': "RdBu",
               'vsymmetric':True,
               #'vmin': -2.0000,
               #'vmax':  2.0000,
               'vmin': -0.5000,
               'vmax':  0.5000,
                },

        'ex': {'title': r"$E_x$",
               'cmap': "PiYG",
               #'vmin': -0.002,
               #'vmax':  0.002,
               #'vmin': -0.0004,
               #'vmax':  0.0004,
               'vsymmetric':True,
                },
        'ey': {'title': r"$E_y$",
               'cmap': "PiYG",
               #'vmin': -0.01,
               #'vmax':  0.01,
               'vsymmetric':True,
                },
        'ez': {'title': r"$E_z$",
               'cmap': "PiYG",
               #'vmin': -0.0001,
               #'vmax':  0.0001,
               #'vmin': -0.000003,
               #'vmax':  0.000003,
               'vsymmetric':True,
                },

        'bx': {'title': r"$B_x$",
               'cmap': "BrBG",
               #'vmin': -0.0001,
               #'vmax':  0.0001,
               #'vmin': -2e-6,
               #'vmax':  2e-6,
               'vsymmetric':True,
                },
        'by': {'title': r"$B_y$",
               'cmap': "BrBG",
               #'vmin': -2e-6,
               #'vmax':  2e-6,
               'vsymmetric':True,
                },
        'bz': {'title': r"$B_z$",
               'cmap': "BrBG",
               'vmin': -0.0005,
               'vmax':  0.0005,
               'vsymmetric':True,
                },
        'dbz': {'title': r"$\delta B_z/B_0$",
               'cmap': "BrBG",
               'vmin': -1.1,
               'vmax':  1.1,
               'derived':True,
               'vsymmetric':True,
                },
        'mpi_grid': {'title': r"mpi",
               'cmap': "tab20",
               #'derived':True,
               'file':'mpi_file',
               },
}


def unpack_args(var):

    #--------------------------------------------------
    # unpack incoming arguments that modify defaults
    args = {}

    #general defaults
    for key in default_values:
        args[key] = default_values[key]

    #overwrite with problem defaults
    try:
        for key in default_problem_values[var]:
            args[key] = default_problem_values[var][key]
    except:
        pass

    #finally, overwrite with user given values
    for key in args:
        try:
            user_val = eval(key)
            if not(user_val == None):
                args[key] = user_val
                print("overwriting {} key with {}".format(key, user_val))
        except:
            pass

    return args

def set_vminmax(val, args, vmin, vmax):

    #print("norm factor: {}".format( norm ))
    print("value at the corner {} / mean val {}".format( val[0,0], np.mean(val)) )
    print("min/max {} / {}".format( np.min(val), np.max(val) ) )
    #print("time of  {}".format(lap2time(info['lap'], conf)))


    #if winsorize
    if not(args['winsorize_min'] == None) or not(args['winsorize_max'] == None):
        wvmin, wvmax = np.quantile(val, [args['winsorize_min'], 1.0-args['winsorize_max']])

        if args['vmin'] == None:
            args['vmin'] = wvmin
        if args['vmax'] == None:
            args['vmax'] = wvmax
    else:
        # else set vmin and vmax using normal min/max
        if args['vmin'] == None:
            args['vmin'] = np.min(val)
        if args['vmax'] == None:
            args['vmax'] = np.max(val)

    # make color limits symmetric
    if args['vsymmetric']:
        vminmax = np.maximum( np.abs(args['vmin']), np.abs(args['vmax']) )
        args['vmin'] = -vminmax
        args['vmax'] =  vminmax


    # finally, re-check that user did not give vmin/vmax
    args['vmin'] = vmin if not(vmin == None) else args['vmin']
    args['vmax'] = vmax if not(vmax == None) else args['vmax']

    #nor the project default
    args['vmin'] = default_problem_values[var]['vmin'] if not(vmin == None) else args['vmin']
    args['vmax'] = default_problem_values[var]['vmax'] if not(vmax == None) else args['vmax']

    return args

 
def plot1d_panel(
        axs, 
        var,
        info, 
        title= None,
        vmin = None,
        vmax = None,
        cmap = None,
        clip = None,
        ):


    #--------------------------------------------------
    # read file

    f5_filename = info['shock_file']

    if not(os.path.isfile(f5_filename)):
        return 0,0,0,0,0

    print('--------------------------------------------------')
    print("reading {}".format(info['shock_file']))

    f5F = h5.File(f5_filename,'r')

    dens = f5F['dens_profile'][()]
    sloc = int( f5F.attrs['shock_loc'][()] )

    ind = np.where(np.isnan(dens))
    dens[ind] = 0
    #ind = np.where(dens > 0)
    #print(dens[ind])

    #--------------------------------------------------
    # box size
    nx, = np.shape(dens)

    walloc = 15.0
    xmin = 0.0 #(walloc)/info['skindepth']
    xmax = (nx-walloc)/info['skindepth']

    xloc = np.arange(len(dens)).astype(float)
    xloc[:] /= info['skindepth']

    #--------------------------------------------------
    # units
    norm = 1.0
    n0 = conf.ppc*2 #*conf.stride**2 #number density per pixel in n_0 
    qe = np.abs(conf.qe)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth
    lenscale = conf.Nx*conf.NxMesh*deltax #(large-eddy size in units of skin depth)

    var = 'rho'
    if var == 'rho' or var == 'logrho':
        #norm = qe*n0 #*conf.gamma
        norm = n0 #*conf.gamma
    dens /= norm

    time = info['lap']*conf.cfl/conf.c_omp

    #--------------------------------------------------
    # find shock front

    sloc2 = np.argmax(dens)
    x_shock_loc = xloc[sloc]

    def where_last(vec):
        a = np.argwhere(vec)
        return a[-1][0] if len(a) else 0

    # v1
    #sloc4 = where_last( dens > conf.shock_density_jump_thr)
    #x_shock_loc = xloc[sloc4]


    # v2
    nu_nd = dens 
    while True:
        sloc4 = where_last( nu_nd > conf.shock_density_jump_thr )
        mean_dens_behind = np.mean(nu_nd[sloc4-3*conf.c_omp:sloc4]) # mean density value behind found location

        #print( 'ind at', ind, 'nu/nd:', nu_nd[ind], 'mean val behind:', mean_dens_behind)

        if mean_dens_behind > conf.shock_density_jump_thr:
            break
        elif sloc4 == 0:
            break # unphysical situation
        else:
            nu_nd = nu_nd[:sloc4] # cut point out and analyze the next

    x_shock_loc = xloc[sloc4]



    #beta_shock
    betash =  x_shock_loc/time # beta_sh
    betashe = 1                # TODO error estimate of beta_sh

    betash = (x_shock_loc-info['x_shock_loc'])/(time-info['prev_time'])


    #print('beta_shock:', betash)

    print('sloc vs argmax', sloc, sloc2, sloc4)
    xloc[:] -= x_shock_loc

    #--------------------------------------------------
    # density jump
    #imin = max(0, sloc4-300) # -200 skind behind front
    #imax = max(imin+30, sloc4-100) # -100skind behind front
    imin = max(0, sloc4-450) # -200 skind behind front
    imax = max(imin+30, sloc4-150) # -100skind behind front
    denslice = dens[imin:imax]

    dens_val = np.mean(denslice)
    dens_err = np.std(denslice)

    #print(denslice)
    print('Nd/Nu:', dens_val, ' +- ', dens_err )




    #--------------------------------------------------
    # plot

    axs[0].plot([time], [x_shock_loc], 'b.')
    #axs[1].plot([time], [betash], 'b.')

    #axs[2].plot([time], [np.mean(denslice)], 'b.')
    axs[2].errorbar([time], [dens_val], yerr=dens_err, color='b')

    return x_shock_loc, betash, betashe, dens_val, dens_err



def plot2d_panel(
        ax, 
        var,
        info, 
        title= None,
        vmin = None,
        vmax = None,
        cmap = None,
        clip = None,
        ):


    print('--------------------------------------------------')
    print("reading {}".format(info['fields_file']))

    args = unpack_args(var)

    f5_filename = info[args['file']]
    f5F = h5.File(f5_filename,'r')

    # normal singular variables
    if not(args['derived']):
        val = read_h5_array(f5F, var)

    # composite quantities
    else:
        print("building composite variable")

        if var == "logrho":
            #val = np.log10(read_h5_array(f5F, 'rho'))
            val = read_h5_array(f5F, 'rho')
        if var == "dbz":
            val = read_h5_array(f5F, 'bz')
            val -= conf.binit
            val /= conf.binit


    # flatten from 3D to 2D
    val = val[:,:,0]


    #--------------------------------------------------
    # get full shape
    nx, ny = np.shape(val)

    #--------------------------------------------------
    # resize to fit aspect ratio
    asp_ratio = nx/ny
    cropr = args['crop_ratio']
    #nxc = int(ny*cropr)
    #print('asp ratio:', asp_ratio)
    #print('new nx {} old nx {} with original aspect ratio {}'.format(nxc, nx, asp_ratio))
    #val = val[:nxc, :] $TODO

    # expand mpi grid array (of tiles) to right size
    if var == "mpi_grid":
        nx *= conf.NxMesh
        ny *= conf.NyMesh

    #walloc = 5.0
    #xmin =-(walloc)/info['skindepth']
    #ymin = 0.0
    #xmax = (nxc-walloc)/info['skindepth']
    #ymax = ny/info['skindepth']

    #x_box = conf.box_shift #- 5.0 #starting point of hires box - 5 cell padding

    xmin = (0.0)/info['skindepth']
    xmax = (nx)/info['skindepth']
    print('xmin and xmax', xmin, xmax)

    ymin = 0.0
    ymax = ny/info['skindepth']

    #--------------------------------------------------
    # normalization
    norm = 1.0
    n0 = conf.ppc*2*conf.stride**2 #number density per pixel in n_0 
 
    qe = np.abs(conf.qe)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth

    lenscale = conf.Nx*conf.NxMesh*deltax #(large-eddy size in units of skin depth)


    if var == 'rho' or var == 'logrho':
        norm = qe*n0 #*conf.gamma
    if var == 'jx' or var == 'jy' or var == 'jz': 
        norm = qe*n0*conf.cfl*conf.cfl
    if var == 'betaperp':
        norm = qe*n0*conf.cfl*conf.cfl
    if var == 'bx' or var == 'by' or var == 'bz': 
        norm = (me_per_qe*conf.cfl**2)/deltax
    if var == 'ex' or var == 'ey' or var == 'ez': 
        norm = (me_per_qe*conf.cfl**2)/deltax
    if var == 'je':
        norm_E = (me_per_qe*conf.cfl**2)/deltax/lenscale
        norm_J = qe*n0*conf.cfl*conf.cfl
        norm = norm_E * norm_J

        #correct for stride size in e/b fields
        norm /= conf.stride**2
        norm /= 1.0e3

    #print('norm: ', norm)
    val = val / norm

    if args['log']:
        val = np.log10(val)

    # set vmin and vmax levels
    args = set_vminmax(val, args, vmin, vmax)


    #--------------------------------------------------
    # print all options
    #print("--- {} ---".format(var))
    #for key in args:
    #    if not(key == None):
    #        print(" setting {}: {}".format(key, args[key]))

    #--------------------------------------------------

    # correct for shock front movement
    xmin -= info['x_shock_loc']
    xmax -= info['x_shock_loc']

    im = imshow(ax, val, xmin, xmax, ymin, ymax,
           cmap = args['cmap'],
           vmin = args['vmin'],
           vmax = args['vmax'],
           clip = args['clip'],
           #aspect=args['aspect'],
           aspect='auto',
           )

    #ax.set_xlim((xmin,xmax))
    #ax.set_ylim((ymin,ymax))

    #xlen = xmax - xmin
    #xs= -20.0 #starting position
    #ax.set_xlim((xs, xs+xlen))

    ax.set_xlim((-20.0, 40.0))
    ax.set_ylim((ymin,ymax))


    #--------------------------------------------------
    # colorbar

    if do_dark:
        plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.97)
    else:
        #ax.set_xlabel(r"$x$ $(c/\omega_p)$")
        #ax.set_ylabel(r"$y$ $(c/\omega_p)$")
        #plt.subplots_adjust(left=0.15, bottom=0.10, right=0.87, top=0.97)
        plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.97)


    wskip = 0.0
    pad = 0.01
    pos = ax.get_position()

    #print(pos)
    axleft   = pos.x0
    axright  = pos.x0 + pos.width
    print('ax pos:', axleft, axright, pos.width, pos.height)

    hpad = 0.08
    axbottom = pos.y0 + hpad*pos.height
    axtop    = pos.y0 + (1.-2.*hpad)*pos.height

    #cax = plt.fig.add_axes([axleft+wskip, axtop+0.01, axright-axleft-2*wskip, 0.01])
    cax = plt.fig.add_axes([axright+pad, axbottom+wskip, 0.01, axtop-axbottom-2*wskip])

    cb = plt.fig.colorbar(im, cax=cax, 
            orientation='vertical',
            ticklocation='right')
    #cb.set_label(args['title']) 

    cax.text(1.0, 1.09, args['title'], transform=cax.transAxes)

    #cb.remove()


#--------------------------------------------------

def plot2d_hires_panel(
        ax, 
        var,
        info, 
        title= None,
        vmin = None,
        vmax = None,
        cmap = None,
        clip = None,
        ):


    print('--------------------------------------------------')
    print("reading {}".format(info['shock_file']))

    args = unpack_args(var)

    #f5_filename = info[args['file']]
    f5_filename = info['shock_file']
    f5F = h5.File(f5_filename,'r')

    # normal singular variables
    if not(args['derived']):
        val = read_h5_array(f5F, var)

    # composite quantities
    else:
        print("building composite variable")

        if var == "logrho":
            #val = np.log10(read_h5_array(f5F, 'rho'))
            val = read_h5_array(f5F, 'rho')
        if var == "dbz":
            val = read_h5_array(f5F, 'bz')
            val -= conf.binit
            val /= conf.binit


    # flatten from 3D to 2D
    val = val[:,:,0]

    #FIXME
    #args['vmin'] = None
    #args['vmax'] = None
    #val[0:240,:] = 0.0


    #val[600:, :] += 0.2

    #--------------------------------------------------
    # get full shape
    nx, ny = np.shape(val)
    print(np.shape(val))

    #--------------------------------------------------
    # resize to fit aspect ratio
    asp_ratio = nx/ny
    cropr = args['crop_ratio']
    nxc = int(ny*cropr)
    #print('asp ratio:', asp_ratio)
    print('new nx {} old nx {} with original aspect ratio {}'.format(nxc, nx, asp_ratio))
    #val = val[:nxc, :] #TODO


    #shift box to start from correct location
    bloc0 = f5F.attrs['bloc0'][()]
    x0 = bloc0/info['skindepth']

    #xmin = x0/info['skindepth']
    #xmax = (nx-bloc0)/info['skindepth']
    #print('xmin xmax', xmin, xmax)
    #print(x0)
    #print(nx,x0)

    #x_box = 400.0 - 5.0 #starting point of hires box - 5 cell padding
    x_box = conf.box_shift #- 5.0 #starting point of hires box - 5 cell padding
    xmin = (0.0+x_box)/info['skindepth']
    xmax = (nx+x_box)/info['skindepth']
    print('xmin and xmax', xmin, xmax)

    ymin = 0.0
    ymax = ny/info['skindepth']


    #--------------------------------------------------
    # normalization
    norm = 1.0
    n0 = conf.ppc*2 #*conf.stride**2 #number density per pixel in n_0 
 
    qe = np.abs(conf.qe)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth

    lenscale = conf.Nx*conf.NxMesh*deltax #(large-eddy size in units of skin depth)


    if var == 'rho' or var == 'logrho':
        norm = qe*n0 #*conf.gamma
    if var == 'jx' or var == 'jy' or var == 'jz': 
        norm = qe*n0*conf.cfl*conf.cfl
    if var == 'betaperp':
        norm = qe*n0*conf.cfl*conf.cfl
    if var == 'bx' or var == 'by' or var == 'bz': 
        norm = (me_per_qe*conf.cfl**2)/deltax
    if var == 'ex' or var == 'ey' or var == 'ez': 
        norm = (me_per_qe*conf.cfl**2)/deltax
    if var == 'je':
        norm_E = (me_per_qe*conf.cfl**2)/deltax/lenscale
        norm_J = qe*n0*conf.cfl*conf.cfl
        norm = norm_E * norm_J

        #correct for stride size in e/b fields
        norm /= conf.stride**2
        norm /= 1.0e3

    #print('norm: ', norm)
    val = val / norm

    if args['log']:
        val = np.log10(val)

    # set vmin and vmax levels
    args = set_vminmax(val, args, vmin, vmax)


    #--------------------------------------------------
    # print all options
    #print("--- {} ---".format(var))
    #for key in args:
    #    if not(key == None):
    #        print(" setting {}: {}".format(key, args[key]))

    #--------------------------------------------------

    # correct for shock front movement
    #xmin -= info['x_shock_loc']
    #xmax -= info['x_shock_loc']

    im = imshow(ax, val, xmin, xmax, ymin, ymax,
           cmap = args['cmap'],
           vmin = args['vmin'],
           vmax = args['vmax'],
           clip = args['clip'],
           #aspect=args['aspect'],
           aspect='auto',
           )

    #ax.set_xlim((xmin,xmax))
    #ax.set_ylim((ymin,ymax))

    #xlen = xmax - xmin
    xlen = 320/info['skindepth']
    xs= -20.0 #starting position

    #ax.set_xlim((xs, xs+xlen))

    ax.set_xlim((-20.0, 40.0))
    ax.set_ylim((ymin,ymax))


    #--------------------------------------------------
    # colorbar

    if do_dark:
        plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.97)
    else:
        #ax.set_xlabel(r"$x$ $(c/\omega_p)$")
        #ax.set_ylabel(r"$y$ $(c/\omega_p)$")
        #plt.subplots_adjust(left=0.15, bottom=0.10, right=0.87, top=0.97)
        plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.97)


    wskip = 0.0
    pad = 0.01
    pos = ax.get_position()

    #print(pos)
    axleft   = pos.x0
    axright  = pos.x0 + pos.width
    print('ax pos:', axleft, axright, pos.width, pos.height)

    hpad = 0.08
    axbottom = pos.y0 + hpad*pos.height
    axtop    = pos.y0 + (1.-2.*hpad)*pos.height

    #cax = plt.fig.add_axes([axleft+wskip, axtop+0.01, axright-axleft-2*wskip, 0.01])
    cax = plt.fig.add_axes([axright+pad, axbottom+wskip, 0.01, axtop-axbottom-2*wskip])

    cb = plt.fig.colorbar(im, cax=cax, 
            orientation='vertical',
            ticklocation='right')
    #cb.set_label(args['title']) 

    cax.text(1.0, 1.09, args['title'], transform=cax.transAxes)

    #cb.remove()


#--------------------------------------------------

    
def quick_build_info(fdir, lap):
    info = {}
    info['lap'] = lap

    info['moms_file']     = fdir + 'moms_'+str(lap)+'.h5'
    info['fields_file']   = fdir + 'flds_'+str(lap)+'.h5'
    info['analysis_file'] = fdir + 'analysis_'+str(lap)+'.h5'
    info['mpi_file']      = fdir + 'mpi_'+str(0)+'.h5'
    info['shock_file']    = fdir + 'shock_'+str(lap)+'.h5'
    info['particle_file'] = fdir + 'test-prtcls'
    
    return info



do_dark = False

if __name__ == "__main__":

    #if do_dark:
    plt.fig = plt.figure(1, figsize=(8,6), dpi=300)
        
    plt.rc('font', family='serif', size=7)
    plt.rc('xtick')
    plt.rc('ytick')

    #else:
    #    plt.fig = plt.figure(1, figsize=(4,3.5), dpi=200)
    #    plt.rc('font',  family='sans-serif')
    #    #plt.rc('text',  usetex=True)
    #    plt.rc('xtick', labelsize=8)
    #    plt.rc('ytick', labelsize=8)
    #    plt.rc('axes',  labelsize=8)
    
    if do_dark:
        plt.style.use('dark_background')

    gs = plt.GridSpec(3, 1)
    gs.update(hspace = 0.5)
    gs.update(wspace = 0.0)
    
    axs = []
    for i in range(3):
        axs.append( plt.subplot(gs[i,0]) )
    
    #--------------------------------------------------
    # command line driven version

    args = pytools.parse_args()
    conf = Configuration_Problem(args.conf_filename)
    fdir = conf.outdir + '/'

    print("plotting {}".format(fdir))
    
    fname_F = "flds"
    fname_A = "analysis"
    fname_P = "test-prtcls"

    #--------------------------------------------------
    # if lap is defined, only plot that one individual round
    Ts_v = [] # time array
    Xs_v = [] # shock front location array

    Rs_v = [] # compression ratio values in time
    Rs_e = [] # compression ratio errors in time

    Vs_v = [] # velocity values in time
    Vs_e = [] # velocity errors in time
    
    offs = 0

    xsh = 0
    prevt = 0
    for lap in range(0, conf.Nt, conf.interval):
        info = quick_build_info(fdir, lap)
        info['skindepth'] = conf.c_omp #/conf.stride

        time = lap*conf.cfl/conf.c_omp
        info['x_shock_loc'] = xsh # store previous shock loc
        info['prev_time'] = prevt  # store previous time

        #print(info['shock_file'])
        if not(os.path.isfile(info['shock_file'])):
            continue

        #--------------------------------------------------
        # analyze shock RH conditions

        xsh, betashv, betashe, comprv, compre = plot1d_panel(axs, 'dens', info)
        prevt = time # update previous time

        #Vs_v.append(betashv) # beta_sh
        #Vs_e.append(betashe) # error of beta_sh
            
        Rs_v.append(comprv) # compression ratio
        Rs_e.append(compre) # error of compression ratio

        Ts_v.append(time)
        Xs_v.append(xsh)

        # estimate betashock slope + error from past points
        if len(Ts_v) > 5:
            xtt = np.array(Ts_v[-5:])
            ytt = np.array(Xs_v[-5:])
            #ind = np.where(ytt > 0)
            #xtt = xtt[ind]
            #ytt = ytt[ind]
            #print(xtt)
            #print(ytt)

            p, V = np.polyfit(xtt, ytt, 1, cov=True)
            #print("x_1: {} +/- {}".format(p[0], np.sqrt(V[0][0])))
            #print("x_2: {} +/- {}".format(p[1], np.sqrt(V[1][1])))

            Vs_v.append( p[0] )
            Vs_e.append( np.sqrt(V[0][0]) ) # error from covariance

            axs[1].errorbar([time], [p[0]], yerr=np.sqrt(V[0][0]), color='b') # plot
            print("betash (5) : {} +/- {}".format(p[0], np.sqrt(V[0][0])))
        else:
            Vs_v.append( betashv)
            Vs_e.append( betashe)

        # longer data sample
        if len(Ts_v) > 20:
            xtt = np.array(Ts_v[-20:])
            ytt = np.array(Xs_v[-20:])
            p, V = np.polyfit(xtt, ytt, 1, cov=True)
            axs[1].errorbar([time], [p[0]], yerr=np.sqrt(V[0][0]), color='r') # plot
            print("betash (20): {} +/- {}".format(p[0], np.sqrt(V[0][0]))) 


    #print('T', Ts_v)
    #print('Rs', Rs_v)
    #print('Re', Rs_e)

    np.savetxt(fdir+'RH.txt', 
               (Ts_v, Rs_v, Rs_e, Vs_v, Vs_e, Xs_v), 
               header='t, R, err, v, err, xloc',
               delimiter=',',
               )


    #axs[0].set_xlim((-400.0, 50.0))
    axs[0].minorticks_on()
    axs[1].minorticks_on()
    axs[2].minorticks_on()

    axs[0].set_xlim((0,10000))
    axs[1].set_xlim((0,10000))
    axs[2].set_xlim((0,10000))

    if True:
        #axs[1].set_ylim((0.2,1.0))
        #axs[2].set_ylim((0.0,5.0))

        axs[1].set_ylim((0.8,1.0))
        axs[2].set_ylim((0.0,5.0))
    else:
        if conf.sigma == 0.0001 and conf.bperp == 1:
            #axs[0].set_ylim((0,1))
            axs[1].set_ylim((0.25,0.4))
            axs[2].set_ylim((3.0,5.0))
        elif conf.sigma == 0.1 and conf.bperp == 1:
            #axs[0].set_ylim((0,1))
            axs[1].set_ylim((0.55,0.6))
            axs[2].set_ylim((3.0,4.0))
        elif conf.sigma == 0.1 and conf.bpar == 1:
            axs[1].set_ylim((0.2,0.4))
            axs[2].set_ylim((3.5,4.5))
    

    axs[0].set_xlabel(r'time $t \omega_p$')
    axs[1].set_xlabel(r'time $t \omega_p$')
    axs[2].set_xlabel(r'time $t \omega_p$')

    axs[0].set_ylabel(r'$x_{sh}$ $(c/\omega_p)$')
    axs[1].set_ylabel(r'$\beta_{sh}$')
    axs[2].set_ylabel(r'$N_d/N_u$')


    fname = fdir +'RH.pdf'

    plt.savefig(fname)

    print('beta fit')

    xtt = np.array(Ts_v[5:])
    ytt = np.array(Xs_v[5:])
    ind = np.where(ytt > 0)
    xtt = xtt[ind]
    ytt = ytt[ind]
    #print(xtt)
    #print(ytt)

    p, V = np.polyfit(xtt, ytt, 1, cov=True)
    print("x_1: {} +/- {}".format(p[0], np.sqrt(V[0][0])))
    print("x_2: {} +/- {}".format(p[1], np.sqrt(V[1][1])))



