import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os
import matplotlib.ticker as ticker
#from matplotlib import cm
import matplotlib 

# 3D
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable

# pytools bindings
import pytools
from init_problem import Configuration_Turbulence as Configuration

from pytools.pybox.box import Box
from pytools.pybox.box import axisEqual3D


default_values = {
    'cmap':"viridis",
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
                'vmax': 4.0,
                'derived':True,
                'file':'moms',
                },
        'jz': {'title': r"$J_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
               'vmin': -1.0000,
               'vmax':  1.0000,
                },
        'bz': {'title': r"$B_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
                },
        'bperp': {'title': r"$B_\perp$",
               'cmap': "magma",
               'vmin':  0.0000,
               'vmax':  10.0000,
               'derived':True,
                },
        'bvec': {'title': r"$B$",
               'cmap': "RdBu",
               'vsymmetric':True,
               'vmin': -1.0000,
               'vmax':  1.0000,
               'derived':True,
                },
        'je': {'title': r"$\mathbf{J} \cdot \mathbf{E}$",
               'cmap': "BrBG",
               'vmin': -0.5,
               'vmax':  0.5,
               'vsymmetric':True,
               'derived':True,
#               'winsorize_min':0.005,
#               'winsorize_max':0.005,
                },
        'gamma': {'title': r"$\Gamma_e$",
                'vmin': 0.0,
                'vmax': 2.0,
                'derived':True,
                'file':'moms',
                },
        'gam_mean': {
            'title': r"$\langle \gamma \rangle$",
            'cmap': "viridis",
            'vmin': 1.0,
            'vmax': 5.0, #rad 2d strong
            #'vmax': 30.0, #rad 2d weak 
            'file':'spec_2d_pixel',
                },
        'gam_top': {
            'title': r"$\log_{10} \gamma_{\mathrm{max}}$",
            'cmap': "inferno",
            'vmin': 0.0,
            'vmax': 1.0,
            #'vmax': 2.5,
            'file':'spec_2d_pixel',
            'log': True,
                },
        'cutoff': {
            'title': r"$\mathcal{A} \gamma_{\mathrm{cool}}$",
            'cmap': "magma",
            'vmin': 0.0,
            #'vmax': 1.0,
            'vmax': 1.2,
            'file':'spec_3d_pixel',
                },
        'mpi': {'title': r"mpi",
               'cmap': "cyclic",
               'derived':True,
               'file':'mpi',
               'vmin': 0,
               #'vmax': 4.0,
               },
}


def get_normalization(var, conf):
    norm = 1.0
    n0 = 2.0*conf.ppc*conf.stride**3 #number density per pixel in n_0 
    qe = np.abs(conf.qe)
    me = np.abs(conf.qe)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth

    lenscale = conf.Nx*conf.NxMesh*deltax/conf.max_mode #(large-eddy size in units of skin depth)

    #print('n0', n0)
    #print('qe', qe)
    #print('me', me)
    #print('stride3', conf.stride**3)

    if var == 'rho': # or var == 'dens':
        norm = n0
    if var == 'dens':
        norm = 2.0*conf.ppc*conf.stride_mom**3  # stide_mom instead of regular stride param
    if var in ['jx','jy','jz']:
        norm = qe*n0*conf.cfl**2
    if var in ["bz", "bx", "by", "bperp"]:
        norm = conf.binit
        #norm = (me_per_qe*conf.cfl**2)/deltax
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

    # bare minimum (visible box sides)
    data = np.zeros((nx,ny,3))
    print('shapes of periphery')
    print(np.shape(data_xy))
    print(np.shape(data_xz))
    print(np.shape(data_yz))

    data[:,:,0] = data_xy[:,:,0]
    data[:,:,1] = data_xz[:,:,0]
    data[:,:,2] = data_yz[:,:,0]

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
    conf = Configuration(args_cli.conf_filename, do_print=False)
    var = args_cli.var
    #lap = args_cli.lap

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

    print(conf.outdir)
    print("plotting {}".format(var))

    #fname_fld = "flds"
    #fname_fld = "moms"
    fname_fld = args['file']
    fname_prtcls = "test-prtcls"

    if args_cli.lap == None:
        laps = np.arange(0, conf.Nt+1, conf.interval)
    else:
        laps = [args_cli.lap]
        

    # read just periphery or full box
    do_periphery = True
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

        print(fname_fld)
        if not(args['derived']):
            data = read_h5(conf.outdir, fname_fld, var, lap)
        else:
            if var == 'gamma':
                vx = read_h5(conf.outdir, fname_fld, "Vxe", lap)
                vy = read_h5(conf.outdir, fname_fld, "Vye", lap)
                vz = read_h5(conf.outdir, fname_fld, "Vze", lap)

                m  = read_h5(conf.outdir, fname_fld, "dense", lap)
                vx[:] /= m
                vy[:] /= m
                vz[:] /= m

                data = 1.0/np.sqrt(1.0 - vx**2 + vy**2 + vz**2)

            if var == 'dens':
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
                data = np.sqrt(bx*bx + by*by) #/conf.binit

            if var == 'mpi':
                data = read_h5(conf.outdir, fname_fld, "mpi_grid", lap)
                args['vmax'] = int(np.max(data))
            #if var == 'densx':
            #    data = read_h5(conf.outdir, fname_fld, "densx", lap)

        if args['log']:
            data = np.log10(data)

        print(np.shape(data))
        nx, ny, nz = np.shape(data)
        #print(data)

        #skim top of due to mom writer bug
        #if args['file'] == 'moms':
        #    data = data[:,:,:-1]

        #--------------------------------------------------
        # normalization
        if True:
            norm = get_normalization(var, conf)
            data = data / norm

        # xxx log density
        #data = np.log10(data)

        print("corner value: {}".format(data[0,0,0]))
        print("min value: {}".format(np.min(data)))
        print("max value: {}".format(np.max(data)))
        print("avg value: {}".format(np.mean(data)))

        # --------------------------------------------------
        # create box
        box = Box(axs[0])

        # box.set_data(np.log10(data))
        box.set_data(data)
        box.vmin = args['vmin']
        box.vmax = args['vmax']
        #cmap = cm.get_cmap(args['cmap'])

        if args['cmap'] == 'cyclic':
            cmap = matplotlib.colormaps['tab20']
            col_list = []
            j = 0
            for i in range(int(args['vmax'])):
                c = cmap(j/20)
                col_list.append(c)
                j += 1
                if j > 20: j = 0
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", col_list)
        else:
            cmap = matplotlib.colormaps[args['cmap']]


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



        axs[0].set_axis_off()
        axs[0].view_init(55.0, -120.0)


    # end of lap loop
    if True:
        # colorbars
        m = plt.cm.ScalarMappable(cmap=cmap)

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
    plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0)

    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    axs[0].set_zlabel("z")

    plt.savefig(fname + ".png")
    #plt.savefig(fname+'.pdf')
