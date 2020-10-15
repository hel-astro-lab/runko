import numpy as np
import h5py as h5
import sys, os

from time import sleep

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

from tvtk.tools import visual

# pytools bindings
import pytools
from problem import Configuration_Packets as Configuration



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
    'vmin':None,
    'vmax':None,
    'vec':False,
    'field':False,
    'isosurf':False,
    'volume':False,
    'volume_slice':False,
}


default_turbulence_values = {
        'bvec': {'title': r"$B$",
               'cmap': "Reds",
               'vmin': 1.0,
               'vmax': 1.4,
               'derived':True,
               'field':True,
               'isosurf':True,
               'vec':True,
                },
        'evec': {'title': r"$E$",
               'cmap': "Reds",
               'vmin': 1.0,
               #'vmax': 1.4,
               'vmax': 3.0,
               'derived':True,
               'field':True,
               'isosurf':True,
               'vec':True,
                },
        'jvec': {'title': r"$J$",
               'cmap': "Reds",
               #'vsymmetric':True,
               'vmin': 0.0,
               'vmax': None,
               'derived':True,
               'field':False,
               'isosurf':True,
               'volume':False,
               'vec':True,
                },
}


def get_normalization(var, conf):
    norm = 1.0
    n0 = 2.0*conf.ppc*conf.stride**3 #number density per pixel in n_0 
    qe = np.abs(conf.qe)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth

    lenscale = conf.Nx*conf.NxMesh*deltax/conf.max_mode #(large-eddy size in units of skin depth)

    if var == 'rho' or var == 'rho2':
        norm = qe*n0
    if var == 'dens':
        norm = n0
    if var == 'jz':
        norm = qe*n0*conf.cfl**2
    if var == 'bz':
        norm = (me_per_qe*conf.cfl**2)/deltax
    if var == 'je':
        norm_E = (me_per_qe*conf.cfl**2)/deltax/lenscale
        norm_J = qe*n0*conf.cfl**2
        norm = norm_E * norm_J

        #correct for stride size in e/b fields
        #norm /= conf.stride**2
        #norm /= 1.0e3

    return norm


def read_h5(outdir, fname_fld, var, lap):
    fields_file = outdir + "/" + fname_fld + "_" + str(lap) + ".h5"
    print(fields_file)
    f5 = h5.File(fields_file, "r")
    return pytools.read_h5_array(f5, var)


if __name__ == "__main__":

    do_print=True
    args_cli = pytools.parse_args()
    conf = Configuration(args_cli.conf_filename, do_print=do_print)
    var = args_cli.var


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
    do_periphery = False
    
    # periphery data is hi-res; take this into account when calculating units
    if do_periphery:
        conf.stride = 1.0


    # plotting 
    #--------------------------------------------------

    fig = mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(600, 600))
    mlab.clf()

    Lx = conf.NxMesh*conf.Nx/conf.stride
    Ly = conf.NyMesh*conf.Ny/conf.stride
    Lz = conf.NzMesh*conf.Nz/conf.stride

    # special branch for pixel spectra files
    if args['file'] == 'spec_3d_pixel':
        print('changing size...')
        Lx = conf.Nx
        Ly = conf.Ny
        Lz = conf.Nz

    box_offset = 0
    for lap in laps:

        if not(args['derived']):
            s = read_h5(conf.outdir, fname_fld, var, lap)
        else:
            if var == 'bvec':
                vx = read_h5(conf.outdir, fname_fld, "bx", lap)/conf.binit
                vy = read_h5(conf.outdir, fname_fld, "by", lap)/conf.binit
                vz = read_h5(conf.outdir, fname_fld, "bz", lap)/conf.binit

            if var == 'evec':
                vx = read_h5(conf.outdir, fname_fld, "ex", lap)/conf.binit
                vy = read_h5(conf.outdir, fname_fld, "ey", lap)/conf.binit
                vz = read_h5(conf.outdir, fname_fld, "ez", lap)/conf.binit

            if var == 'jvec':
                vx = read_h5(conf.outdir, fname_fld, "jx", lap)
                vy = read_h5(conf.outdir, fname_fld, "jy", lap)
                vz = read_h5(conf.outdir, fname_fld, "jz", lap)


        if args['vec']:
            #3D vector field
            s = np.sqrt(vx**2 + vy**2 + vz**2)
            src = mlab.pipeline.vector_field(vx, vy, vz)
            magnitude = mlab.pipeline.extract_vector_norm(src)
        else:
            #scalar
            src = mlab.pipeline.scalar_field(s)

        #print scalar min max
        print("min/max:", np.min(s), np.max(s))

        # get min/max if not given
        if args['vmin'] == None:
            args['vmin'] = np.min(s)
        if args['vmax'] == None:
            args['vmax'] = np.max(s)
        
        # 3D vector field
        if args['field']:

            #extract vector objects
            src = mlab.pipeline.vector_field(vx, vy, vz)
            magnitude = mlab.pipeline.extract_vector_norm(src)

            # streamline tubes
            streamline = mlab.pipeline.streamline(
                    magnitude, 
                    seedtype='sphere',
                    integration_direction='both',
                    colormap='YlGnBu', #default colormap for field line intensity
                    vmin=args['vmin'], 
                    vmax=args['vmax'],
                    #opacity=0.8,
                    seed_resolution=10,
                    seed_scale=1.0,
                    seed_visible=False,
                    )

            streamline.stream_tracer.maximum_propagation = 500.0

            streamline.streamline_type = 'tube'
            streamline.tube_filter.vary_radius = True
            streamline.tube_filter.number_of_sides = 5
            streamline.tube_filter.radius = 0.2
            streamline.tube_filter.radius_factor = 5.0

            #mlab.pipeline.vectors(
            #        src, 
            #        mask_points=30, 
            #        scale_factor=5.,
            #        opacity=0.8,
            #        colormap='plasma',
            #        )


        ############################################################ 
        # iso surfaces
        if args['isosurf']:

            dv = args['vmax']-args['vmin']
            contours=np.array([
                #0.025,
                0.125,
                0.25,
                0.5,
                0.75,
                ])
            contours = args['vmin'] + dv*contours

            surf = mlab.pipeline.iso_surface(
                    magnitude, 
                    contours=contours.tolist(),
                    opacity=0.2,
                    colormap=args['cmap'],
                    transparent=True,
                    )


        ############################################################ 
        # volume rendering
        if args['volume']:

            vol = mlab.pipeline.volume(
                    src,
                    vmin=args['vmin'], 
                    vmax=args['vmax'],
                    #colormap=args['cmap'],
                    )

        if args['volume_slice']:
            mlab.pipeline.image_plane_widget(
                    src,
                    plane_orientation='x_axes',
                    slice_index=Lx//2,
                    vmin=args['vmin'], 
                    vmax=args['vmax'],
                    colormap=args['cmap'],
                    opacity=0.8,
                    )

            mlab.pipeline.image_plane_widget(
                    src,
                    plane_orientation='y_axes',
                    slice_index=Ly//2,
                    vmin=args['vmin'], 
                    vmax=args['vmax'],
                    colormap=args['cmap'],
                    opacity=0.8,
                    )

            mlab.pipeline.image_plane_widget(
                    src,
                    plane_orientation='z_axes',
                    slice_index=Lz//4,
                    vmin=args['vmin'], 
                    vmax=args['vmax'],
                    colormap=args['cmap'],
                    opacity=0.8,
                    )


        ############################################################ 
        #draw axis and simulation box outline

        extent=[0, Lx, 0, Ly, 0, Lz,]
        mlab.outline(extent=extent, opacity=0.8,)

        Lxp = conf.NxMesh*conf.Nx/conf.c_omp
        Lyp = conf.NyMesh*conf.Ny/conf.c_omp
        Lzp = conf.NzMesh*conf.Nz/conf.c_omp

        axes = mlab.axes(
                extent=extent,
                ranges=(
                    -Lxp//2, Lxp//2, 
                    -Lyp//2, Lyp//2, 
                     0, Lz), 
                xlabel=r'x',
                ylabel=r'y',
                zlabel=r'z',
                opacity=0.8,
                nb_labels=5
                )

        axes.label_text_property.font_family = 'times'
        axes.label_text_property.font_size = 4

        axes.title_text_property.font_family = 'times'
        axes.title_text_property.font_size = 7

        ax=mlab.axes()
        print("format:", ax.axes.label_format)
        ax.axes.label_format='%3.0f'


        # normalize units
        #norm = get_normalization(var, conf)
        #data = data / norm
    
        if False:
            s = data**2
            src = mlab.pipeline.scalar_field(s)

            vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(data), 
                    vmin=1.5, vmax=3.0 )#, color='YlGn')


        #vol.module_manager.scalar_lut_manager.lut_mode = 'YlGnBu'
        #mlab.pipeline.iso_surface(src, contours=[s.min()+0.2*s.ptp(), ], opacity=0.1)
        #mlab.pipeline.iso_surface(src, contours=[3.0, 4.0, ], opacity=0.05)
        #mlab.pipeline.iso_surface(src, contours=[s.max()-0.1*s.ptp(), ],)
        #mlab.pipeline.image_plane_widget(src,
        #                            plane_orientation='z_axes',
        #                            slice_index=10,
        #                        )

    #--------------------------------------------------
    
    mlab.view(50.0, 70.0, 6.25*Lx)
    #mlab.view(50.0, 70.0, 3.25*Lx)
    #mlab.view(50.0, 70.0, 250)
    mlab.draw()

    slap = str(lap).rjust(4, "0")
    fname = conf.outdir + "/3d_render_" + var + "_" + slap
    mlab.savefig(fname+".png")
    
    mlab.show()
