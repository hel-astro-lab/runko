# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os
import matplotlib.ticker as ticker

from init_problem import Configuration_Shocks as Configuration_Problem

#from combine_files import get_file_list
#from combine_files import combine_tiles

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
        ax, 
        var,
        info, 
        offs = 0,
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

    if not(os.path.isfile(f5_filename)):
        return 0

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

    walloc = 5.0
    xmin = 0.0 #(walloc)/info['skindepth']
    xmax = (nx-walloc)/info['skindepth']

    xloc = np.arange(len(dens)).astype(float)
    xloc[:] /= info['skindepth']

    if sloc >= len(xloc):
        sloc = -1

    x_shock_loc = xloc[sloc]
    xloc[:] -= x_shock_loc

    xlen = 50.0
    xs= -20.0 #starting position
    #ax.set_xlim((xs, xs+xlen))
    #ax.set_xlim((xmin,xmax))


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

    #--------------------------------------------------
    # plot

    dens[:] += offs
    ax.plot(xloc, dens, 'b-', lw=0.1)

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
    plt.fig = plt.figure(1, figsize=(8,8), dpi=300)
        
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

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.5)
    gs.update(wspace = 0.0)
    
    axs = []
    for i in range(1):
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
    
    offs = 0
    for lap in range(0, conf.Nt, conf.interval):
        info = quick_build_info(fdir, lap)
        info['skindepth'] = conf.c_omp #/conf.stride
        print(info['shock_file'])
        
        offs += 1.0

        #--------------------------------------------------
        # plot shock density profile and get shock front position
        xsh = plot1d_panel(axs[0], 'dens', info, offs)
        info['x_shock_loc'] = xsh


    #TODO: use this to set xlim
    axs[0].set_xlim((-120.0, 120.0))

    #axs[0].set_ylim((0.0, 4.0))
    axs[0].minorticks_on()

    slap = str(args.lap).rjust(7, '0')
    fname = fdir +'dens_mnt.pdf'
    plt.savefig(fname)

    #wide
    axs[0].set_xlim((-500.0, 500.0))
    fname = fdir +'dens_mnt_wide.pdf'
    plt.savefig(fname)

