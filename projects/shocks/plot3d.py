import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os
import matplotlib.ticker as ticker
from matplotlib import cm

# 3D
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable

# pytools bindings
import pytools
from problem import Configuration_Problem as Configuration

from pytools.pybox.box import Box
from pytools.pybox.box import axisEqual3D



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
    'file':'flds',
    'log':False,
    'vmin':-1,
    'vmax':+1,
}


default_turbulence_values = {
        'rho': {'title': r"$n/n_0$",
                'vmin': 0.0,
                'vmax': 4.0,
                },
        'dens': {'title': r"$n/n_0$",
                'vmin': 0.0,
                'vmax': 6.0,
                'derived':True,
                'file':'moms',
                },
        'logdens': {'title': r"$\log_{10} n/n_0$",
                'vmin':-2.0,
                'vmax': 1.0,
                'derived':True,
                'file':'moms',
                'log':'True',
                },
        'jz': {'title': r"$J_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
               'vmin': -2.0000,
               'vmax':  2.0000,
                },
        'jperp': {'title': r"$J_\perp$",
               'cmap': "magma",
               #'vsymmetric':True,
               'vmin': 0.0000,
               'vmax': 0.1000,
                'derived':True,
                },
        'bdens': {'title': r"$U_B$",
               'cmap': "magma",
               #'vsymmetric':True,
               'vmin': 0.0000,
               'vmax': 1.0000,
                'derived':True,
                },
        'bz': {'title': r"$B_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
                },
        'bperp': {'title': r"$B_\perp$",
               'cmap': "magma",
               'vmin': -1.0000,
               'vmax':  1.0000,
               'derived':True,
                },
        'gamma': {'title': r"$\Gamma_e$",
                'vmin': 0.0,
                'vmax': 20.0,
                'derived':True,
                'file':'moms',
                },
}



def get_normalization(var, conf):
    norm = 1.0
    n0 = 2.0*conf.ppc*conf.stride**3 #number density per pixel in n_0 
    qe = np.abs(conf.qe)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth

    #lenscale = conf.Nx*conf.NxMesh*deltax/conf.max_mode #(large-eddy size in units of skin depth)
    lenscale = 1.0 #one skin depth

    if var == 'rho': # or var == 'dens':
        norm = qe*n0
    if var == 'dens' or var == 'logdens':
        norm = n0
    if var == 'jz':
        norm = qe*n0*conf.cfl**2
    if var == 'bz':
        norm = (me_per_qe*conf.cfl**2)/deltax
    if var == 'je':
        norm_E = (me_per_qe*conf.cfl**2)/deltax/lenscale
        norm_J = qe*n0*conf.cfl**2
        norm = norm_E * norm_J

    return norm


def read_full_box(outdir, fname_fld, var, lap):
    fields_file = outdir + "/" + fname_fld + "_" + str(lap) + ".h5"
    f5 = h5.File(fields_file, "r")
    return pytools.read_h5_array(f5, var)


def read_periphery(outdir, fname, var, lap):

    #fields_file = conf.outdir + "/" + fname_fld + "_" + str(lap) + ".h5"
    fname = "slices" #overwrite
    xy_file = outdir + "/" + fname + "-xy_" + str(lap) + ".h5"
    xz_file = outdir + "/" + fname + "-xz_" + str(lap) + ".h5"
    yz_file = outdir + "/" + fname + "-yz_" + str(lap) + ".h5"

    f5_xy = h5.File(xy_file, "r")
    f5_xz = h5.File(xz_file, "r")
    f5_yz = h5.File(yz_file, "r")

    data_xy = pytools.read_h5_array(f5_xy, var)
    data_xz = pytools.read_h5_array(f5_xz, var)
    data_yz = pytools.read_h5_array(f5_yz, var)

    nx,ny,ns0 = np.shape(data_xy)
    print(nx,ny,ns0)
    nx,nz,ns1 = np.shape(data_xz)
    print(nx,nz,ns1)
    ny,nz,ns2 = np.shape(data_yz)
    print(ny,nz,ns2)

    # full data
    #data = np.zeros((nx,ny,nz))

    #xy
    #data[:,:,0 ] = data_xy[:,:,0]
    #data[:,:,-1] = data_xy[:,:,0]

    ##xz
    #data[:,0, :] = data_xz[:,:,0]
    #data[:,-1,:] = data_xz[:,:,0]

    ##yz
    #data[0, :,:] = data_yz[:,:,0]
    #data[-1,:,:] = data_yz[:,:,0]

    # bare minimum (visible box sides)
    data = np.zeros((nx,ny,3))
    data[:,:,0] = data_xy[:,:,0]
    data[:,:,1] = data_xz[:,:,0]
    #data[:,:,2] = data_yz[:,:,0]

    return data


def load_args(var):

    #general defaults
    args = {}
    for key in default_values:
        args[key] = default_values[key]

    #overwrite with turbulence defaults
    try:
        for key in default_turbulence_values[var]:
            args[key] = default_turbulence_values[var][key]
    except:
        pass

    return args

def load_data(var, args, lap):
    fname_fld = args['file']

    if not(args['derived']):
        data = read_h5(conf.outdir, fname_fld, var, lap)
    else:
        if var == 'gamma':
            vx = read_h5(conf.outdir, fname_fld, "Vxe", lap)
            vy = read_h5(conf.outdir, fname_fld, "Vye", lap)
            vz = read_h5(conf.outdir, fname_fld, "Vze", lap)
            data = 1.0/np.sqrt(1.0 - vx**2 + vy**2 + vz**2)
        if var == 'dens' or var == 'logdens':
            de = read_h5(conf.outdir, fname_fld, "dense", lap)
            dp = read_h5(conf.outdir, fname_fld, "densp", lap)
            data = de + dp
        if var == 'je':
            jx = read_h5(conf.outdir, fname_fld, "jx", lap)
            jy = read_h5(conf.outdir, fname_fld, "jy", lap)
            jz = read_h5(conf.outdir, fname_fld, "jz", lap)
            ex = read_h5(conf.outdir, fname_fld, "ex", lap)
            ey = read_h5(conf.outdir, fname_fld, "ey", lap)
            ez = read_h5(conf.outdir, fname_fld, "ez", lap)
            data = jx*ex + jy*ey + jz*ez
        if var == 'bperp':
            bx = read_h5(conf.outdir, fname_fld, "bx", lap)
            by = read_h5(conf.outdir, fname_fld, "by", lap)
            data = bx*bx + by*by
        if var == 'jperp':
            jy = read_h5(conf.outdir, fname_fld, "jy", lap)
            jz = read_h5(conf.outdir, fname_fld, "jz", lap)
            data = np.sqrt(jy*jy + jz*jz)
        if var == 'bdens':
            bx = read_h5(conf.outdir, fname_fld, "bx", lap)
            by = read_h5(conf.outdir, fname_fld, "by", lap)
            bz = read_h5(conf.outdir, fname_fld, "bz", lap)
            data = bx*bx + by*by + bz*bz

    return data

if __name__ == "__main__":

    do_print = True
    do_dark = False

    if do_dark:
        fig = plt.figure(1, figsize=(8, 4.0), dpi=300)

        plt.rc("font", family="serif", size=7)
        plt.rc("xtick")
        plt.rc("ytick")

        plt.style.use("dark_background")
    else:
        fig = plt.figure(1, figsize=(3.45, 3.45), dpi=200)
        plt.rc("font", family="sans-serif")
        # plt.rc('text',  usetex=True)
        plt.rc("xtick", labelsize=8)
        plt.rc("ytick", labelsize=8)
        plt.rc("axes", labelsize=8)
        # plt.style.use('dark_background')

    gs = plt.GridSpec(1, 1)
    gs.update(hspace=0.0)
    gs.update(wspace=0.0)

    axs = []
    # axs.append( plt.subplot(gs[0,0]) )
    axs.append(fig.add_subplot(111, projection="3d"))

    # --------------------------------------------------
    # prepare data

    args_cli = pytools.parse_args()
    conf = Configuration(args_cli.conf_filename, do_print=do_print)
    var = args_cli.var

    args = load_args(var)

    print(conf.outdir)
    print("plotting {}".format(var))

    #fname_fld = args['file']
    #fname_prtcls = "test-prtcls"

    if args_cli.lap == None:
        laps = np.arange(0, conf.Nt+1, conf.interval)
    else:
        laps = [args_cli.lap]

    # read just periphery or full box
    do_periphery = False
    if not(do_periphery):
        def read_h5(outdir, fname, var, lap):
            return read_full_box(outdir, fname, var, lap)
    else:
        def read_h5(outdir, fname, var, lap):
            return read_periphery(outdir, fname, var, lap)
    
    # periphery data is hi-res; take this into account when calculating units
    if do_periphery:
        conf.stride = 1.0


    box_offset = 0
    for lap in laps:
        data = load_data(var, args, lap)

        # limit box size
        data = data[1:,:,:] #cut out reflector

        #xlen = 800.0
        #xlen = 1300.0
        xlen = 700.0
        xleni = np.int(xlen*conf.c_omp/conf.stride) #into units of cells
        data = data[0:xleni,:,:]

        #manual striding
        stride = 5
        #data = data[::stride,::stride,::stride]


        print(np.shape(data))
        nx, ny, nz = np.shape(data)

        #XXX TODO remove
        #skim top of due to mom error
        #if args['file'] == 'moms':
        #    data = data[:,:,:-1]

        #--------------------------------------------------
        # normalization
        if True:
            norm = get_normalization(var, conf)
            data = data / norm

        if args['log']:
            data = np.log10(data)

        print("corner value: {}".format(data[0,0,0]))
        print("min value: {}".format(np.min(data)))
        print("max value: {}".format(np.max(data)))

        # --------------------------------------------------
        # create box
        box = Box(axs[0])

        # box.set_data(np.log10(data))
        box.set_data(data)
        box.vmin = args['vmin']
        box.vmax = args['vmax']
        cmap = cm.get_cmap(args['cmap'])

        box_offset = 0.0
        box.z0 -= box_offset

        aspect_ratio = nx/ny
        print("box aspect ratio: {}".format(aspect_ratio))
        box.dx = aspect_ratio
        box.dy = 1.2
        box.dz = 1.0 #1.0


        # surface rendering
        if do_periphery:
            print("drawing only periphery using reduced data")
            box.draw_top(cmap=cmap,   data_slice=data[:,:,0])
            box.draw_front(cmap=cmap, data_slice=data[:,:,1])
            box.draw_left(cmap=cmap,  data_slice=data[:,:,2])
            box.draw_outline()
        else:
            box.draw_left(cmap=cmap)
            box.draw_front(cmap=cmap)
            box.draw_top(cmap=cmap)
            box.draw_outline()



        # draw ticks
        if False:
            #ticks = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
            ticks = [100, 200, 300, 400]

            L = 426.0 #box length
            tick_locations = np.array([0, 100, 200, 300, 400])
            print(tick_locations)

            # top x axis
            #box.set_ticklabels(
            #        11,
            #        ticks,
            #        position='top',
            #        along='x',
            #        direction='out',
            #        )

            # bottom x axis
            box.set_ticklabels(
                    L,
                    tick_locations,
                    [],
                    position='top',
                    along='x',
                    direction='in',
                    )

            # bottom x axis
            box.set_ticklabels(
                    L,
                    tick_locations,
                    ticks,
                    label='$x$ ($c/\omega_p$)',
                    position='bottom',
                    along='x',
                    direction='in',
                    offs_label = 0.18,
                    rotation=20.0,
                    )


            # bottom y axis
            box.set_ticklabels(
                    L,
                    tick_locations,
                    [],
                    position='top',
                    along='y',
                    direction='in',
                    )
            # bottom y axis
            box.set_ticklabels(
                    L,
                    tick_locations,
                    ticks,
                    label='$y$ ($c/\omega_p$)',
                    position='bottom',
                    along='y',
                    direction='in',
                    offs_label = 0.20,
                    rotation=-20.0,
                    )

            # bottom z axis
            box.set_ticklabels(
                    L,
                    tick_locations,
                    ticks,
                    label='$z$ ($c/\omega_p$)',
                    position='top',
                    along='z',
                    direction='in',
                    offs_label = 0.22,
                    rotation=90.0,
                    )

            # bottom z axis
            box.set_ticklabels(
                    L,
                    tick_locations,
                    [],
                    position='bottom',
                    along='z',
                    direction='in',
                    )


        # slices; same data
        if True:
            for loc in np.linspace(0.05, 0.75, 8):
                locx = loc*aspect_ratio
                print("slize at:", locx)
                box.draw_exploded_slice("left-front", locx, 3.9, cmap=cmap)


        #back exploded panels
        if True:
            #var2 = 'rho'
            var2 = 'bdens'
            print("plotting: ", var2)

            args = load_args(var2)
            data = load_data(var2, args, lap)

            data = data[6:,:,:] #cut off reflector
            data = data[0:xleni,:,:] #limit box length

            #stride down
            #data = data[::stride,::stride,::stride]

            box.set_data(data)

            box.vmin = args['vmin']
            box.vmax = args['vmax']
            cmap = cm.get_cmap(args['cmap'])

            off_bot    = 1.7
            off_back   = 1.7
            off_left   = 0.7
            off_right  = 0.7

            box.draw_exploded_panels_outline("bottom", off=off_bot)
            box.draw_exploded_panels_outline("back",   off=off_back)
            
            box.draw_exploded_bottom(off=off_bot,   cmap=cmap)
            box.draw_exploded_back( off=off_back,   cmap=cmap)

            #box.draw_exploded_panels_outline("left",   off=off_left)
            #box.draw_exploded_left(  off=off_left,  cmap=cmap)

            #box.draw_exploded_panels_outline("right",  off=off_right)
            #box.draw_exploded_right( off=off_right, cmap=cmap)


            #off = -1.95 #this puts them in front of the box
            #box.draw_exploded_panels_outline("bottom", off=off)
            #box.draw_exploded_panels_outline("left",   off=off)
            #box.draw_exploded_panels_outline("right",  off=off)
            #
            #box.draw_exploded_bottom(off=off, cmap=cmap)
            #box.draw_exploded_left(  off=off, cmap=cmap)
            #box.draw_exploded_right( off=off, cmap=cmap)




        axs[0].set_axis_off()
        # axs[0].view_init(45.0, -110.0)

        #axs[0].view_init(55.0, -120.0)
        axs[0].view_init(35.0, -140.0)


    # end of lap loop

    if False:
        # colorbars
        m = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap(args['cmap']))

        m.set_array( [args['vmin'], args['vmax'] ])
        m.set_clim(vmin=args['vmin'], vmax=args['vmax'])
        cbaxes = fig.add_axes([0.3, 0.06, 0.4, 0.01])  # [left, bottom, width, height],
        cb = plt.colorbar(
            m, cax=cbaxes, 
            orientation="horizontal", 
            ticklocation="bottom",
        )
        cb.ax.tick_params(labelsize=6) 

        for axis in ['top','bottom','left','right']:
          cb.ax.spines[axis].set_linewidth(0.5)

        # cb.set_label(r'$J$', rotation=0)
        fig.text(0.22, 0.06, args['title'], fontsize=8, ha='center', va='center')


    # force square box
    axisEqual3D(axs[0])

    # axs[0].set_title('Step {}'.format(lap+1))

    slap = str(lap).rjust(4, "0")
    fname = conf.outdir + "/3d_" + var + "_" + slap
    #plt.subplots_adjust(left=-0.45, bottom=-0.45, right=1.35, top=1.45)
    plt.subplots_adjust(left=-0.45, bottom=-0.7, right=1.45, top=1.6)

    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    axs[0].set_zlabel("z")

    plt.savefig(fname + ".png")
    #plt.savefig(fname+'.pdf')
