import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats
from scipy.optimize import curve_fit

from visualize import imshow

from configSetup import Configuration
from combine_files import get_file_list
from combine_files import combine_tiles

from scipy.ndimage.filters import gaussian_filter

from parser import parse_input
import argparse



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
    'aspect':1,
    'vsymmetric':None,
    'winsorize_min':0.005,
    'winsorize_max':0.005,
    'title':'',
    'derived':False,
}

default_shock_values = {
        'rho': {'title': r"$\rho$",
                'vmin': 0.0,
                #'vmax': 0.8,
                },
        'jz': {'title': r"$J_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
               #'vmin': -2.0000,
               #'vmax':  2.0000,
                },
        'bz': {'title': r"$B_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
                },
}



def read_var(f5F, var):
    try:
        val = f5F[var][:,:,0]
    except:
        nx = f5F['Nx'].value
        ny = f5F['Ny'].value
        nz = f5F['Nz'].value
        #print("reshaping 1D array into multiD with {} {} {}".format(nx,ny,nz))

        val = f5F[var][:]
        val = np.reshape(val, (ny, nx, nz))
        val = np.transpose(val[:,:,0])

    return val


def build_bdens(f5F):
    bx = read_var(f5F, "bx")
    by = read_var(f5F, "by")
    bz = read_var(f5F, "bz")

    return bx*bx + by*by + bz*bz

def build_edens(f5F):
    ex = read_var(f5F, "ex")
    ey = read_var(f5F, "ey")
    ez = read_var(f5F, "ez")

    return ex*ex + ey*ey + ez*ez


def nandiv(x, y):
    return x/y



def plot2d_shock_single(
        ax, 
        var,
        info, 
        title= None,
        vmin = None,
        vmax = None,
        cmap = None,
        clip = None,
        ):

    #--------------------------------------------------
    # unpack incoming arguments that modify defaults
    args = {}

    #general defaults
    for key in default_values:
        args[key] = default_values[key]

    #overwrite with shock defaults
    try:
        for key in default_shock_values[var]:
            args[key] = default_shock_values[var][key]
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

    print('--------------------------------------------------')
    print("reading {}".format(info['fields_file']))

    f5F = h5.File(info['fields_file'],'r')

    # normal singular variables
    if not(args['derived']):
        val = read_var(f5F, var)

    # composite quantities
    else:
        print("building composite variable")

        if var == "bdens":
            val = build_bdens(f5F)
        if var == "edens":
            val = build_edens(f5F)

    #--------------------------------------------------
    # get shape
    nx, ny = np.shape(val)
    #print("nx={} ny={}".format(nx, ny))

    walloc = 5.0
    xmin =-(walloc)/info['skindepth']
    ymin = 0.0
    xmax = (nx-walloc)/info['skindepth']
    ymax = ny/info['skindepth']


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
    args['vmin'] = default_shock_values[var]['vmin'] if not(vmin == None) else args['vmin']
    args['vmax'] = default_shock_values[var]['vmax'] if not(vmax == None) else args['vmax']

    #--------------------------------------------------
    # print all options
    print("--- {} ---".format(var))
    for key in args:
        if not(key == None):
            print(" setting {}: {}".format(key, args[key]))

    #--------------------------------------------------

    im = imshow(ax, val, xmin, xmax, ymin, ymax,
           cmap = args['cmap'],
           vmin = args['vmin'],
           vmax = args['vmax'],
           clip = args['clip'],
           aspect=args['aspect'],
           )

    #cax, cb = colorbar(im)

    #--------------------------------------------------
    # colorbar
    #ax.set_xlim((200.0, 600.0))
    ax.set_xlim((0.0, 400.0))
    ax.set_ylim((0.0, 64.0))


    if do_dark:
        plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.97)
    else:
        ax.set_xlabel(r"$x$ $(c/\omega_p)$")
        ax.set_ylabel(r"$y$ $(c/\omega_p)$")
        plt.subplots_adjust(left=0.15, bottom=0.10, right=0.87, top=0.97)

    #axleft    = 0.10
    #axbottom  = 0.06
    #axright   = 0.96
    #axtop     = 0.92
    wskip = 0.0
    pad = 0.01
    pos = ax.get_position()
    #print(pos)
    axleft   = pos.x0
    axbottom = pos.y0
    axright  = pos.x0 + pos.width
    axtop    = pos.y0 + pos.height

    print(axleft)
    print(axbottom)
    print(axright)
    print(axtop)

    #cax = plt.fig.add_axes([axleft+wskip, axtop+0.01, axright-axleft-2*wskip, 0.01])
    cax = plt.fig.add_axes([axright+pad, axbottom+wskip, 0.01, axtop-axbottom-2*wskip])

    cb = plt.fig.colorbar(im, cax=cax, 
            orientation='vertical',
            ticklocation='right')
    #cb.set_label(args['title']) 

    cax.text(1.0, 1.03, args['title'], transform=cax.transAxes)


    #if do_dark:
    #    for axis in ['top', 'bottom', 'left', 'right']:
    #        cax.spines[axis].set_visible(False)
    #        #cax.spines[axis].set_linewidth(0)
    #        #cax.spines[axis].set_color('red')
    #ax.set_title(title)

    #plot2d_particles(ax, info, args)
    if do_prtcls:
        plot2d_particles(ax, info, args)

        # if plotting with particles; add a tag to name
        varp = var + '_p'
    else:
        varp = var

    slap = str(info['lap']).rjust(4, '0')
    if do_dark:
        fname = fdir + varp+'_{}.png'.format(slap)
        plt.savefig(fname)
    else:
        fname = fdir + varp+'_{}.pdf'.format(slap)
        plt.savefig(fname)
    cb.remove()



class TestParticles:

    def __init__(self):
        self.prtcls = {}
        self.times = []

    def add_file(self, info):
        fname = info['particle_file'] + "_" + str(info['lap']) + ".h5"
        try:
            #print("adding ", fname)
            f5F = h5.File(fname,'r')
        except:
            #print("failed")
            return

        nx = f5F['Nx'].value
        ny = f5F['Ny'].value
        nz = f5F['Nz'].value
        #print("reading particles")
        #print("reshaping 1D array into multiD with ({} {} {}) {}".format(nx,ny,nz, info['lap']))

        xloc = np.reshape( f5F['x'][:]    ,(nx,-1))
        yloc = np.reshape( f5F['y'][:]    ,(nx,-1)) 
        zloc = np.reshape( f5F['z'][:]    ,(nx,-1)) 
        ux   = np.reshape( f5F['vx'][:]   ,(nx,-1)) 
        uy   = np.reshape( f5F['vy'][:]   ,(nx,-1)) 
        uz   = np.reshape( f5F['vz'][:]   ,(nx,-1)) 
        wgt  = np.reshape( f5F['wgt'][:]  ,(nx,-1)) 
        ids  = np.reshape( f5F['id'][:]   ,(nx,-1)) 
        procs= np.reshape( f5F['proc'][:] ,(nx,-1))

        #print("non zero prtcls:", np.count_nonzero(xloc))

        nt, npp = np.shape(xloc)
        for i in range(npp):
            key = (ids[0,i], procs[0,i])

            if not(key in self.prtcls):
                self.prtcls[key] = []

            for t in range(nt):
                self.prtcls[key].append([
                        xloc[t,i],
                        yloc[t,i],
                        zloc[t,i],
                        ux[t,i],
                        uy[t,i],
                        uz[t,i],
                        wgt[t,i]])

        #print("adding {}".format(info['lap']))
        self.times.append(info['lap'])



# Plot test particle locations in a 2D plane
def plot2d_particles(ax, info, args):

    print('--------------------------------------------------')
    print("reading {}".format(info['fields_file']))

    f = open(prtcl_file, 'r')
    xx = f.read().splitlines()
    f.close()

    marker_style = dict(
            color='white', 
            linestyle=':', 
            marker='o',
            markersize=3, 
            mew=1,
            fillstyle='none',
            markerfacecoloralt='black',
            alpha=0.8,
            )

    marker_history_style = dict(
            color='white', 
            linestyle='solid', 
            linewidth=1.0,
            alpha=0.5,
            )

    stride=2.0
    #np = 100
    #xx = prtcls.xloc[:,::1000]/stride
    #yy = prtcls.yloc[:,::1000]/stride
    #ax.plot(xx, yy, **marker_style)

    for line in xx:
        k1,k2, maxgam = line.split(',')
        key = (int(k1),int(k2))
        #print("getting {}".format(key))

        xs = []
        ys = []
        zs = []
        #vx = []
        #vy = []
        #vz = []
        gs = []

        for prtcl in prtcls.prtcls[key]:
            #print(prtcl)
            xs.append( prtcl[0]/stride )
            ys.append( prtcl[1]/stride )
            zs.append( prtcl[2]/stride )

            ux = prtcl[3]
            uy = prtcl[4]
            uz = prtcl[5]
            gam = np.sqrt(1.0 + ux*ux + uy*uy + uz*uz)
            gs.append( gam )

        tsteps = 5
        #ax.plot(xs[-tsteps:], ys[-tsteps:], **marker_history_style)

        #plot head
        marker_style['markersize']=np.log10(gs[-1])*3.0
        marker_style['mew']=np.log10(gs[-1])*2.0
        ax.plot(xs[-1], ys[-1], **marker_style)



        #plot history by wrapping around periodic boundaries
        xarr = [ xs.pop() ]
        yarr = [ ys.pop() ]
        garr = [ gs.pop() ]
        ds = 0.0

        tstep = 5
        ic = 0
        while ic < tstep:
            if len(xs) > 0:
                ds = np.sqrt( (xarr[-1] - xs[-1])**2 + (yarr[-1] - ys[-1])**2 )
                if ds > 100:
                    ax.plot( xarr, yarr, **marker_history_style)
                    
                    xarr = [ xs.pop() ]
                    yarr = [ ys.pop() ]
                    garr = [ gs.pop() ]
                    ic += 1
                else:
                    xarr.append( xs.pop() )
                    yarr.append( ys.pop() )
                    garr.append( gs.pop() )
                    ic += 1
            else:
                break

        ax.plot( xarr, yarr, **marker_history_style)




#--------------------------------------------------

def build_info(fdir, lap):
    info = {}
    info['lap'] = lap
    info['fields_file']   = fdir + 'fields_'+str(args.lap)+'.h5'
    info['analysis_file'] = fdir + 'analysis'+str(args.lap)+'.h5'
    
    return info
    
def quick_build_info(fdir, lap):
    info = {}
    info['lap'] = lap
    info['fields_file']   = fdir + 'flds_'+str(args.lap)+'.h5'
    info['analysis_file'] = fdir + 'analysis_'+str(args.lap)+'.h5'
    info['particle_file'] = fdir + 'test-prtcls'
    
    return info



# global particle container
prtcls = TestParticles()

do_dark = True

if __name__ == "__main__":

    if do_dark:
        plt.fig = plt.figure(1, figsize=(8,2.0), dpi=300)
        
        plt.rc('font', family='serif', size=7)
        plt.rc('xtick')
        plt.rc('ytick')
    else:
        plt.fig = plt.figure(1, figsize=(4,2.5), dpi=200)
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
    
    #ax.set_xlabel(r"$x$ $(c/\omega_p)$")
    #ax.set_ylabel(r"$y$ $(c/\omega_p)$")

    #--------------------------------------------------
    # command line driven version

    conf, fdir, args = parse_input()

    fdir += '/'
    print("plotting {}".format(fdir))
    
    fname_F = "flds"
    fname_A = "analysis"
    fname_P = "test-prtcls"

    if do_dark:
        plt.style.use('dark_background')

    do_prtcls = False
    prtcl_file = fdir + "sorted_top10_prtcls.txt"


    # if lap is defined, only plot that one individual round
    if not(args.lap == None):

        #info = build_info(fdir, args.lap)
        info = quick_build_info(fdir, args.lap)
        info['skindepth'] = conf.c_omp/conf.stride
        if do_prtcls:
            for hlaps in range(args.lap-5*conf.interval, args.lap+1, conf.interval):
                #print("adding lap from history {}".format(hlaps)) 
                info['lap'] = hlaps
                prtcls.add_file(info)

        plot2d_shock_single(axs[0], args.var, info)

    # else plot every file that there is
    # FIXME
    else:
    
        files_F = get_file_list(fdir, fname_F)
        files_A = get_file_list(fdir, fname_A)
        files_P = get_file_list(fdir, fname_P)
    
        for lap, f in enumerate(files_F):

            #info = build_info(fdir, lap)
            info = {}
            info['lap'] = lap*conf.interval
            info['fields_file'  ]   = files_F[lap]
            info['skindepth'] = conf.c_omp/conf.stride
            #info['particle_file'  ] = files_P[lap]

            if do_prtcls:
                info['particle_file'] = fdir + 'test-prtcls'
                prtcls.add_file(info)

            try:
                info['analysis_file'] = files_A[lap]
            except:
                print("no analysis file found...")
                pass

            plot2d_shock_single(axs[0], args.var, info)

    
