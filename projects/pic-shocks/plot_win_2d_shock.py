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
                'vmax': 5.0,
                },
        'logrho': {'title': r"$\log_{10} n/n_0$",
                'vmin': -1.0,
                'vmax':  1.0,
                'derived':True,
                'log':True,
                },
        'jx': {'title': r"$J_x$",
               'cmap': "RdBu",
               'vsymmetric':True,
               #'vmin': -2.0000,
               #'vmax':  2.0000,
               #'vmin': -0.5000,
               #'vmax':  0.5000,
               #'vmin': -0.0001,
               #'vmax':  0.0001,
                },
        'jy': {'title': r"$J_y$",
               'cmap': "RdBu",
               'vsymmetric':True,
               #'vmin': -2.0000,
               #'vmax':  2.0000,
               #'vmin': -0.5000,
               #'vmax':  0.5000,
               #'vmin': -0.0001,
               #'vmax':  0.0001,
                },
        'jz': {'title': r"$J_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
               #'vmin': -0.0001,
               #'vmax':  0.0001,
               #'vmin': -0.500,
               #'vmax':  0.500,
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
               #'vmin': -0.0005,
               #'vmax':  0.0005,
               'vsymmetric':True,
                },
        'dbz': {'title': r"$\delta B_z/B_0$",
               'cmap': "BrBG",
               #'vmin': -4.0,
               #'vmax':  4.0,
               'derived':True,
               'vsymmetric':True,
                },
        'mpi_grid': {'title': r"mpi",
               'cmap': "tab20",
               #'derived':True,
               'file':'mpi_file',
               },
        'gam_0': {'title': r"$\mathrm{log}_{10}(\gamma_e-1)$",       },
        'ux_0': {'title': r"$\beta_{x,e}\gamma_e$", },
        'uy_0': {'title': r"$\beta_{y,e}\gamma_e$", },
        'uz_0': {'title': r"$\beta_{z,e}\gamma_e$", },
        'gam_1':{'title': r"$\gamma_i - 1$",        },
        'ux_1': {'title': r"$\beta_{x,p}\gamma_p$", },
        'uy_1': {'title': r"$\beta_{y,p}\gamma_p$", },
        'uz_1': {'title': r"$\beta_{z,p}\gamma_p$", },
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
                #print("overwriting {} key with {}".format(key, user_val))
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

    #--------------------------------------------------
    # read file

    f5_filename = info['shock_file']
    f5F = h5.File(f5_filename,'r')

    #dens = read_h5_array(f5F, 'dens_profile')
    dens = f5F['dens_profile'][()]
    #print(f5F['dens_profile'].keys())
    #sloc = f5F[0].attr['shock_loc']

    sloc = int( f5F.attrs['shock_loc'][()] )


    #--------------------------------------------------
    # box size
    #dens = dens[:512]
    nx, = np.shape(dens)
    print('nx: ',nx)

    walloc = 15.0
    xmin = 0.0 #(walloc)/info['skindepth']
    xmax = (nx-walloc)/info['skindepth1d']

    xloc = np.arange(len(dens)).astype(float)
    xloc[:] /= info['skindepth1d']

    #XXX

    x_shock_loc = xloc[sloc]
    xloc[:] -= x_shock_loc

    xlen = 50.0
    xs= -20.0 #starting position
    #ax.set_xlim((xs, xs+xlen))
    #ax.set_xlim((xmin,xmax))
    #ax.set_xlim((-20.0, 40.0))
    #ax.set_xlim((-600.0, 200.0))

    ax.set_ylim((0.0, 5.0))
    ax.minorticks_on()

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
        norm = n0 #*conf.gamma
    dens /= norm

    #--------------------------------------------------
    # plot
    ax.plot(xloc, dens, 'b-')

    #ax.axvline(x=x_shock_loc)
    print('sloc', sloc)
    print(xloc[sloc])

    #pos = ax.get_position()
    ##print(pos)
    ##axleft   = pos.x0
    ##axright  = pos.x0 + pos.width
    #axleft   = 0.193333333
    #axright  = 0.806666666
    #axwidth = axright-axleft
    #axheight = 0.076666666
    #axbottom = pos.y0
    #ax.set_position([axleft, axbottom, axwidth, axheight])
    #print('ax pos:', axleft, axright)

    return x_shock_loc



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
            if conf.use_maxwell_split:
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
        norm = n0 #*conf.gamma
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

    im = pytools.visualize.imshow(ax, val, xmin, xmax, ymin, ymax,
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

    ax.set_xlim((-600.0, 200.0))
    ax.set_ylim((ymin,ymax))


    #--------------------------------------------------
    # colorbar

    if do_dark:
        plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.97)
    else:
        ax.set_xlabel(r"$x$ $(c/\omega_p)$")
        ax.set_ylabel(r"$y$ $(c/\omega_p)$")
        plt.subplots_adjust(left=0.15, bottom=0.10, right=0.87, top=0.97)


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

    #shift box to start from correct location
    bloc0 = f5F.attrs['bloc0'][()]
    x0 = bloc0/info['skindepth']

    x_box = conf.box_shift/conf.box_stride #- 5.0 #starting point of hires box - 5 cell padding
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
        norm = n0 #*conf.gamma
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

    print('mean before', np.mean(val[np.nonzero(val)]))
    print('norm: ', norm)
    val = val / norm
    print('mean after', np.mean(val[np.nonzero(val)]))

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

    im = pytools.visualize.imshow(ax, val, xmin, xmax, ymin, ymax,
           cmap = args['cmap'],
           vmin = args['vmin'],
           vmax = args['vmax'],
           clip = args['clip'],
           #aspect=args['aspect'],
           aspect='auto',
           )

    #ax.set_xlim((xmin,xmax))
    #ax.set_ylim((ymin,ymax))

    #ax.set_xlim((-600.0, 200.0))
    ax.set_ylim((ymin,ymax))


    #--------------------------------------------------
    # colorbar
    if do_dark:
        plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.97)
    else:
        ax.set_xlabel(r"$x$ $(c/\omega_p)$")
        ax.set_ylabel(r"$y$ $(c/\omega_p)$")
        plt.subplots_adjust(left=0.15, bottom=0.10, right=0.87, top=0.97)

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

def plot2d_pspec(
        ax, 
        var,
        info, 
        title= None,
        vmin = None,
        vmax = None,
        cmap = None,
        clip = None,
        extent = None,
        ):

    print('--------------------------------------------------')
    print("reading {}".format(info['shock_file']))

    args = unpack_args(var)

    f5_filename = info['shock_file']
    f5F = h5.File(f5_filename,'r')
    val = f5F[var][()]

    if var == 'gam_0' or var == 'gam_1':
        bins = f5F['bin_gams'][()]
    else:
        bins = f5F['bin_betagam'][()]

    sloc = int( f5F.attrs['shock_loc'][()] )
    nx, ny = np.shape(val)

    #log of histogram
    val = np.log10(val)
    
    #shift box to start from correct location
    bloc0 = f5F.attrs['bloc0'][()]
    x0 = bloc0/info['skindepth1d']

    x_box = conf.box_shift #starting point of hires box 

    xmin = (0.0+x_box)/info['skindepth1d']
    xmax = (nx+x_box)/info['skindepth1d']

    #manually set args
    args['cmap'] = 'plasma_r'
    args['vmin'] =-2.0
    args['vmax'] = 5.0

    if var == 'gam_0' or var == 'gam_1':
        ymin = np.log10(bins[0])
        ymax = np.log10(bins[-1])
    else:
        ymin = bins[0]
        ymax = bins[-1]

    args = set_vminmax(val, args, vmin, vmax)
    print('vmin/max:', args['vmin'], args['vmax'])

    im = pytools.visualize.imshow(ax, val, xmin, xmax, ymin, ymax,
           cmap = args['cmap'],
           vmin = args['vmin'],
           vmax = args['vmax'],
           clip = None, #args['clip'],
           #aspect=args['aspect'],
           aspect='auto', #XXX
           )

    ax.set_ylim((ymin, ymax))
    print('ymin/max:', ymin, ymax)

    uminmax = 40
    if var == 'ux_0' or var == 'ux_1':
        ax.set_ylim((-uminmax, uminmax))
    if var == 'uy_0' or var == 'uy_1':
        ax.set_ylim((-uminmax, uminmax))
    if var == 'uz_0' or var == 'uz_1':
        ax.set_ylim((-uminmax, uminmax))

    ax.minorticks_on()

    # set xlim from previous panel
    #axleft1, axbottom1, axright1, axtop1, xmin, xmax = extent

    #ax.set_xlim((-600, 200))

    if not(var == 'gam_0' or var == 'gam_1'):
        ax.axhline(y=0.0, linestyle='dotted', alpha=0.6, color='w')

    #if var == 'gam_0' or var == 'gam_1':
    #    ax.set_yscale('log')

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
    #cax.text(1.0, 1.05, args['title'], transform=cax.transAxes)

    #cax.text(1.0, 1.09, args['title'], transform=cax.transAxes)
    #cax.text(1.05, 1.05, args['title'], transform=cax.transAxes)
    ax.set_xlabel(r"$x$ $(c/\omega_p)$")
    ax.set_ylabel(args['title'], color='w') 

    return


def quick_build_info(fdir, lap):
    info = {}
    info['lap'] = lap

    info['moms_file']     = fdir + 'moms_'+str(args.lap)+'.h5'
    info['fields_file']   = fdir + 'flds_'+str(args.lap)+'.h5'
    info['analysis_file'] = fdir + 'analysis_'+str(args.lap)+'.h5'
    info['mpi_file']      = fdir + 'mpi_'+str(0)+'.h5'
    info['shock_file']    = fdir + 'shock_'+str(args.lap)+'.h5'
    info['particle_file'] = fdir + 'test-prtcls'
    
    return info



do_dark = True

if __name__ == "__main__":

    if do_dark:
        plt.fig = plt.figure(1, figsize=(8,8), dpi=300)
        
        plt.rc('font', family='serif', size=7)
        plt.rc('xtick', top = True)
        plt.rc('ytick', right = True)
    
    if do_dark:
        plt.style.use('dark_background')

    gs = plt.GridSpec(16, 1)
    gs.update(hspace = 0.1)
    gs.update(wspace = 0.0)
    
    axs = []
    for i in range(16):
        axs.append( plt.subplot(gs[i,0]) )
    
    #--------------------------------------------------
    # command line driven version

    args = pytools.parse_args()
    conf = Configuration_Problem(args.conf_filename)

    if args.lap == None:
        print('Specify lap with --lap xxxx')

    fdir = conf.outdir + '/'
    print("plotting {}".format(fdir))
    
    fname_F = "flds"
    fname_A = "analysis"
    fname_P = "test-prtcls"

    #--------------------------------------------------
    # if lap is defined, only plot that one individual round

    info = quick_build_info(fdir, args.lap)
    info['skindepth']   = conf.c_omp/conf.box_stride  #/conf.stride
    info['skindepth1d'] = conf.c_omp                  #/conf.stride

    #--------------------------------------------------
    # plot shock density profile and get shock front position
    xsh = plot1d_panel(axs[0], 'dens', info)
    info['x_shock_loc'] = xsh

    plot2d_hires_panel(axs[1], 'rho',  info)

    plot2d_hires_panel(axs[2], 'ex', info)
    plot2d_hires_panel(axs[3], 'ey', info)
    plot2d_hires_panel(axs[4], 'ez', info)

    plot2d_hires_panel(axs[5], 'bx', info)
    plot2d_hires_panel(axs[6], 'by', info)
    plot2d_hires_panel(axs[7], 'dbz', info)

    plot2d_hires_panel(axs[8], 'jx', info)
    plot2d_hires_panel(axs[9], 'jy', info)
    plot2d_hires_panel(axs[10],'jz', info)

    # NOTE: mpi grid is always from first lap
    # this is done in already during construction of 'info'
    #plot2d_hires_panel(axs[11], 'jz',  info)
    plot2d_panel(axs[11], 'mpi_grid', info)

    plot2d_pspec(axs[12],  'gam_0', info, ) #extent=ext)
    plot2d_pspec(axs[13],  'ux_0',  info, ) #extent=ext)
    plot2d_pspec(axs[14],  'uy_0',  info, ) #extent=ext)
    plot2d_pspec(axs[15], 'uz_0',   info, ) #extent=ext)


    for ax in axs:
        ax.set_xlim((-50.0, 50.0))
        #ax.set_xlim((-200.0, 200.0))

    slap = str(args.lap).rjust(7, '0')
    if do_dark:
        fname = fdir +'movwin_shock_{}.png'.format(slap)
    else:
        fname = fdir +'movwin_shock_{}.pdf'.format(slap)

    plt.savefig(fname)


