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
from init_problem import Configuration_Shocks as Configuration

from pytools.pybox.box import Box
from pytools.pybox.box import axisEqual3D

import pyvista as pv

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
    'file':'shock', #'flds',
    'log':False,
    'vmin':-1,
    'vmax':+1,
}


default_turbulence_values = {
        'rho': {'title': r"$n/n_0$",
                'vmin': 0.0,
                'vmax': 4.0,
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
               'vmin': -1.0000,
               'vmax':  1.0000,
               'derived':True,
                },
        'bvec': {'title': r"$B$",
               'cmap': "RdBu",
               'vsymmetric':True,
               'vmin': -1.0000,
               'vmax':  1.0000,
               'derived':True,
                },
}


def get_normalization(var, conf):

    norm = 1.0
    n0 = 2.0*conf.ppc #*conf.stride**3 #number density per pixel in n_0 
    qe = np.abs(conf.qe)
    me = np.abs(conf.me)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth

    lenscale = conf.Nx*conf.NxMesh*deltax #/conf.max_mode #(large-eddy size in units of skin depth)

    if var == 'rho': # or var == 'dens':
        norm = 2.0*conf.ppc #*conf.gamma #*me
    if var == 'dens':
        norm = n0
    if var in ['jx', 'jy', 'jz']:
        norm = qe*n0*conf.cfl**2
    if var in ['bx', 'by', 'bz']:
        norm = conf.binit
        #norm = (me_per_qe*conf.cfl**2)/deltax
    if var == 'je':
        norm_E = (me_per_qe*conf.cfl**2)/deltax/lenscale
        norm_J = qe*n0*conf.cfl**2
        norm = norm_E * norm_J

        #correct for stride size in e/b fields
        #norm /= conf.stride**2
        norm /= 1.0e3

    return norm


def read_full_box(outdir, fname_fld, var, lap, stride=1):
    fields_file = outdir + "/" + fname_fld + "_" + str(lap) + ".h5"
    f5 = h5.File(fields_file, "r")
    return pytools.read_h5_array(f5, var, stride=stride)


do_print = True

if __name__ == "__main__":

    args_cli = pytools.parse_args()
    conf = Configuration(args_cli.conf_filename, do_print=do_print)
    var = args_cli.var # NOTE uncomment for basic functionality 
    #var = 'bvec' # manually set the plotted variable

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


    fname_fld = args['file']
    fname_prtcls = "test-prtcls"


    lap = args_cli.lap
    #lap = 2000

    # reading funciton for data
    def read_h5(outdir, fname, var, lap, stride=1):
        return read_full_box(outdir, fname, var, lap, stride=stride)

    #--------------------------------------------------
    print(fname_fld)
    s = 1
    rho= read_h5(conf.outdir, fname_fld, "rho", lap,stride=s).T

    jx = read_h5(conf.outdir, fname_fld, "jx", lap, stride=s).T
    jy = read_h5(conf.outdir, fname_fld, "jy", lap, stride=s).T
    jz = read_h5(conf.outdir, fname_fld, "jz", lap, stride=s).T

    bx = read_h5(conf.outdir, fname_fld, "bx", lap, stride=s).T
    by = read_h5(conf.outdir, fname_fld, "by", lap, stride=s).T
    bz = read_h5(conf.outdir, fname_fld, "bz", lap, stride=s).T

    # normalize
    rho/= get_normalization('rho',conf)

    jx /= get_normalization('jx', conf)
    jy /= get_normalization('jy', conf)
    jz /= get_normalization('jz', conf)

    bx /= get_normalization('bx', conf)
    by /= get_normalization('by', conf)
    bz /= get_normalization('bz', conf)


    dx = conf.stride/conf.c_omp # skindepth resolution
    origin = 0,0,0

    print(np.shape(rho))
    nx, ny, nz = np.shape(rho)

    Lx = nx*dx
    Ly = ny*dx
    Lz = nz*dx

    #--------------------------------------------------
    # zoom
    #   z  y  x
    if True:

        # point 1
        #w = 100
        #xcen = 0.55*Lx
        #ycen = 0.38*Ly
        #midz = int(Lz*0.4)

        #TODO shift location a bit 

        # point2
        #w = 100
        #xcen = 0.40*Lx
        #ycen = 0.55*Ly

        #ymin = int(ycen - w//2)
        #ymax = int(ycen + w//2)

        #xmin = int(xcen - w//2)
        #xmax = int(xcen + w//2)

        #jx[ :, ymin:ymax, xmin:xmax] = 10
        #jy[ :, ymin:ymax, xmin:xmax] = 10 
        #jz[ :, ymin:ymax, xmin:xmax] = 10 

        xmin = 500
        xmax = 1500

        rho = rho[:, :, xmin:xmax]
        jx  = jx[ :, :, xmin:xmax]
        jy  = jy[ :, :, xmin:xmax]
        jz  = jz[ :, :, xmin:xmax]
        bx  = bx[ :, :, xmin:xmax]
        by  = by[ :, :, xmin:xmax]
        bz  = bz[ :, :, xmin:xmax]

    print(np.shape(rho))
    nz, ny, nx = np.shape(rho)

    #-------------------------------------------------- 
    p = pv.Plotter(
            #lighting='three lights'
            off_screen=True,
            )

    Lx = nx*dx
    Ly = ny*dx
    Lz = nz*dx

    midx = Lx//2
    midy = Ly//2
    midz = Lz//2

    Lx_box = 1.0*Lx
    Ly_box = 1.0*Ly


    #--------------------------------------------------
    #if True: # density
    if var == 'dens': # volume rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # scalar field
        clim = (0.0, 4.0)

        print('min/max', np.min(rho), np.max(rho))

        #ind = np.where(rho < 1.5)
        #rho[ind] = 0.0

        mesh['rho'] = np.ravel(rho)
        p.add_mesh( mesh.outline_corners(), color="k" )

        ops = np.zeros(256)
        #ops[:] = 128*np.linspace(0.0, 1.0, 256)**2.0
        ops[:] = 128*np.linspace(0.0, 1.0, 256) #**2.0

        #ops[0:192] = 0.0
        #ops[192:]  = np.linspace(0.0, 128, 64)
        #ops[64:]  = np.linspace(0.0, 256, 192)

        p.add_volume(mesh,
                     scalars='rho',
                     opacity=ops, #'geom',
                     clim=clim,
                     show_scalar_bar=False,
                     cmap='viridis',
                     shade=False,
                    )


    #--------------------------------------------------
    # b field
    #if True: # streamlines
    if var == 'dens': # streamlines

        m1 = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        vectors = np.empty((m1.n_points, 3))
        vectors[:,0] = np.ravel(bx)
        vectors[:,1] = np.ravel(by)
        vectors[:,2] = np.ravel(bz)

        m1['vectors'] = vectors

        print(m1)

        #p.add_mesh(m1.outline(), color="k")

        #mesh = pv.ImageData(
        #    dimensions=(nx, ny, nz), 
        #    spacing=(dx, dx, dx), 
        #    origin=origin)
        b2 = np.sqrt( bx**2 + by**2 + bz**2 )
        print('min/max', np.min(b2), np.max(b2), np.mean(b2))
        m1['scalar'] = np.ravel(1.0 + b2)

        aspect_rat = nx/ny
        print('aspect rat', aspect_rat)
    
        seed = pv.Plane(center=(midx,midy,midz), 
                        i_size=0.98*Lx, 
                        j_size=0.98*Ly,
                        i_resolution=int(10*aspect_rat),
                        j_resolution=10,
                        direction=(0,0,1),
                        )

        #p.add_mesh(seed.outline(), color='r')

        stream = m1.streamlines_from_source(
            seed,
            vectors='vectors',
            max_time=1e5,
            #initial_step_length=0.1,
            integration_direction='both',
            #integration_direction='forward',
            )


        #stream, src = m1.streamlines(
        #    'vectors', 
        #    return_source=True, 
        #    terminal_speed=0.0, 
        #    n_points=200, 
        #    source_radius=200,
        #    source_center=(0, 0, 0),
        #)

        p.add_mesh(
                stream.tube(
                    radius=0.025,
                    scalars="scalar",
                    radius_factor=100.0,
                    n_sides=30,
                    ),
                show_scalar_bar=False,
                cmap='reds',
                #cmap='inferno',
                ambient=0.1,
                clim=(0.0, 3.0),
                   )

        #ribbon
        #p.add_mesh(
        #        stream.ribbon(
        #            width=0.05,
        #            scalars="scalar",
        #            factor=300.0,
        #            ),
        #        show_scalar_bar=False,
        #        cmap='reds',
        #        #cmap='inferno',
        #        ambient=0.1,
        #        clim=(1.0, 3.0),
        #           )


    #--------------------------------------------------
    #if False: # volume rendering
    if var == 'cur': # volume rendering

        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # scalar field
        clim = (-0.3, 0.3)
        
        dbz = bz - 1.0 # -B_0
        dbz[:,:, :520] = 0.0
        
        print('min/max dbz', np.min(dbz), np.max(dbz))
        mesh['dbz'] = np.ravel(dbz)

        ops = np.zeros(256)
        ops[0:128] = np.linspace(128, 0, 128)
        ops[128:]  = np.linspace(0, 128, 128)

        p.add_volume(mesh,
                     scalars='dbz',
                     opacity=ops, #'linear', #ops, #'geom',
                     clim=clim,
                     show_scalar_bar=False,
                     cmap='BrBG',
                     shade=False,
                    )


    #--------------------------------------------------
    #if False: # try isocountours rendering
    if var == 'cur': # try isocountours rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # j scalar field
        #j = np.sign(jz)*np.sqrt( jx**2 + jy**2 + jz**2 )
        j = jy

        clim = (-2.0, 2.0)
        mesh['j'] = np.ravel(j)
        p.add_mesh( mesh.outline_corners(), color="k" )

        ops = np.zeros(256)
        ops[0:128] = np.linspace(128, 0, 128)
        ops[128:]  = np.linspace(0, 128, 128)

        #ops = np.zeros(256)
        #ops[0:64]     = np.linspace(128.0, 0.0,   64)
        #ops[64:128]   = 0.0
        #ops[128:192:] = 0.0
        #ops[192:]     = np.linspace(0.0,   128.0, 64)

        print(ops)

        p.add_volume(mesh,
                     scalars='j',
                     opacity=ops, #'linear',
                     clim=clim,
                     show_scalar_bar=False,
                     #cmap='seismic',
                     cmap='RdBu',
                     shade=False,
                    )


    #--------------------------------------------------
    #--------------------------------------------------
    #--------------------------------------------------


    if False: # current slice
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # j scalar field
        #jz = np.sign(jz)*np.sqrt( jx**2 + jy**2 + jz**2 )
        mesh['j'] = np.ravel(jz)

        #mesh = mesh.gaussian_smooth(std_dev=0.7)

        cmap = 'RdBu'
        clim = (-1.5, 1.5)


        # y slice
        if True: 
            single_slice = mesh.slice(
                    normal=[0,1,0],
                    origin=(midx,midy,midz),
                    )

            p.add_mesh(single_slice, 
                       cmap=cmap,
                       clim=clim,
                       show_scalar_bar=False,
                       opacity=1.0,
                       )

        # z slice
        if True: 
            single_slice2 = mesh.slice(
                    normal=[0,0,1],
                    origin=(midx,midy,midz),
                    )

            p.add_mesh(single_slice2, 
                       cmap=cmap,
                       clim=clim,
                       show_scalar_bar=False,
                       opacity=1.0,
                       )




    #--------------------------------------------------
    # zoom in b field
    if False: # streamlines

        m1 = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        vectors = np.empty((m1.n_points, 3))
        vectors[:,0] = np.ravel(bx)
        vectors[:,1] = np.ravel(by)
        vectors[:,2] = np.ravel(bz) #*0.8

        m1['vectors'] = vectors

        print(m1)

        #p.add_mesh(m1.outline(), color="k")
        

        seed = pv.Plane(center=(midx,midy,midz), 
                        i_size=Lx_box, 
                        j_size=Ly_box,
                        i_resolution=30,
                        j_resolution=6,
                        direction=(0,0,1),
                        )

        #p.add_mesh(seed.outline(), color='r')

        b2 = np.sqrt( bx**2 + by**2 + bz**2 )
        print('min/max', np.min(b2), np.max(b2), np.mean(b2))
        m1['scalar'] = np.ravel(1.0 + b2)

        #stream = m1.streamlines_from_source(
        #    seed,
        #    vectors='vectors',
        #    max_time=2e4,
        #    initial_step_length=0.01,
        #    integration_direction='both',
        #    )

        stream = m1.streamlines(
             'vectors',
             return_source=False, 
        #    terminal_speed=0.0, 
            n_points=150, 
            source_radius=40,
            source_center=(midx, midy, midz),
        )

        p.add_mesh(
                stream.tube(
                    #radius=0.4,
                    scalars="scalar",
                    radius=0.04,
                    radius_factor=18.0,
                    ),
                show_scalar_bar=False,
                cmap='viridis',
                ambient=0.5,
                opacity=1.0,
                clim=(0.5, 3.0),
                   )




    if False: # rho contours rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # rho
        print('min/max', np.min(rho), np.max(rho), ' mean', np.mean(rho))

        mesh['scalar'] = np.ravel(rho)
        mesh = mesh.gaussian_smooth(std_dev=0.3)

        p.add_mesh( mesh.outline_corners(), color="k")

        contours = mesh.contour(np.linspace(1.5, 4.0, 10))
        #contours = contours.smooth_taubin(50)
        p.add_mesh(contours,
                    opacity=0.05,
                    clim=(0.0, 10.0),
                    show_scalar_bar=False,
                    )

        hc  = mesh.contour([4.0])
        ##hc = hc.smooth_taubin(10, pass_band=10.0)
        p.add_mesh(hc, 
                   opacity=1.0,
                   clim=(0.,10.0),
                   show_scalar_bar=False,
                   )



    if False: # try isocountours rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # j scalar field
        clim = (0, 10.0)

        ind = np.where(rho < 2.0)
        rho[ind] = 0.0

        mesh['rho'] = np.ravel(rho)
        p.add_mesh( mesh.outline_corners(), color="k" )


        ops = np.zeros(256)
        #ops[0:192] = 0.0
        #ops[192:]  = np.linspace(0.0, 128, 64)

        ops[64:]  = np.linspace(0.0, 256, 192)

        p.add_volume(mesh,
                     scalars='rho',
                     opacity='linear',
                     clim=clim,
                     show_scalar_bar=False,
                     cmap='viridis',
                     shade=False,
                    )



    if False: # try isocountours rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # j scalar field
        jz = np.sign(jz)*np.sqrt( jx**2 + jy**2 + jz**2 )
        mesh['scalar'] = np.ravel(jz)
        #mesh = mesh.gaussian_smooth(std_dev=0.8)

        p.add_mesh( mesh.outline_corners(), color="k")

        vv = 0.8
        clim = (-1.0, 1.0)

        contours = mesh.contour([0.5*vv, 0.7*vv, 0.9*vv])
        p.add_mesh(contours,
                    opacity=0.3,
                    clim=clim,
                    show_scalar_bar=False,
                    cmap='RdBu',
                    )

        contours = mesh.contour([-0.5*vv, -0.7*vv, -0.9*vv])
        p.add_mesh(contours,
                    opacity=0.3,
                    clim=clim,
                    show_scalar_bar=False,
                    cmap='RdBu',
                    )

        hcp = mesh.contour([+vv])
        p.add_mesh(hcp, 
                   opacity=1.0,
                    clim=clim,
                    show_scalar_bar=False,
                    cmap='RdBu',
                   )

        hcm = mesh.contour([-vv])
        p.add_mesh(hcm, 
                   opacity=1.0,
                    clim=clim,
                   show_scalar_bar=False,
                   cmap='RdBu',
                   )


    if False: # try isocountours rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # j scalar field
        jabs = np.sqrt(jx**2 + jy**2 + jz**2)
        mesh['scalar'] = np.ravel(jabs)

        #mesh = mesh.gaussian_smooth(std_dev=0.8)

        p.add_mesh( mesh.outline_corners(), color="k")

        contours = mesh.contour(np.linspace(0.2, 1.0, 5))
        #contours = contours.smooth_taubin(50)
        p.add_mesh(contours,
                    opacity=0.1,
                    clim=(0.3,2.0),
                    show_scalar_bar=False,
                    )

        hc  = mesh.contour([1.0])
        #hc = hc.smooth_taubin(10, pass_band=10.0)

        p.add_mesh(hc, 
                   opacity=1.0,
                   clim=(0.3,2.0),
                   show_scalar_bar=False,
                   )


        #contours = mesh.contour([1.0], method='marching_cubes')
        ##co = grid.contour([0], values, method='flying_edges')
        #dist = np.linalg.norm(mesh.points, axis=1)

        #mesh.plot(
        #        scalars=dist, 
        #        smooth_shading=True, 
        #        specular=1, 
        #        cmap="plasma", 
        #        show_scalar_bar=False)

    #cpos = [(-670.0201105669008, 595.4204616543602, 918.8903921685189),
    #     (520.0390507995008, 8.315342200304187, 17.7006156267516),
    #     (0.36242940747691593, 0.9238260923231535, -0.12324883666332916)]

    #cpos = [(-515.313099088259, 931.4812834769558, 823.3663162550966),
    #     (525.9195970918121, 4.151508661196132, 30.28989917383742),
    #     (0.4263045917217798, 0.8147566848345564, -0.39298338335650984)]
    cpos = [(-494.7043995313236, -850.17820536772, 929.309505034927),
        (479.8400327799021, 69.62919237730318, 47.666866155907435),
        (0.40148406453978863, 0.3753699092261755, 0.8354088682604126)]
    p.camera_position = cpos

    #p.enable_depth_peeling(100)
    #p.enable_anti_aliasing('msaa')
    #p.enable_ssao(kernel_size=32)

    #old
    #p.save_graphic("3d_jz_b_v2.pdf", title="", raster=True, painter=True)

    cpos = p.show(return_cpos=True, auto_close=False)



    slap = str(lap).rjust(5, '0')
    p.screenshot(conf.outdir + "/" + "3d_" + var + "_" + slap + ".png", scale=2)

    print('camera pos:', cpos)

    print('azimuth', p.camera.azimuth)
    print('elevation', p.camera.elevation)
    print('view_angle', p.camera.view_angle)
    print('up', p.camera.up)

    if False:
        p.open_movie(conf.outdir + "/" + "3d_jz_b.mp4", quality=8)

        for az in np.linspace(0, 360, 360):
            p.camera.azimuth = az
            #p.screenshot("3d_jz_b_{}.png".format(int(az)) )
            p.write_frame()

        #for el in np.linspace(0, 50, 100):
        #    p.camera.elevation = el
        #    p.write_frame()

        ## zoom in
        #for va in np.linspace(30, 10, 100):
        #    p.camera.view_angle = va
        #    p.write_frame()

        #for pause in range(20):
        #    p.write_frame()

        ## zoom out
        #for va in np.linspace(10, 30, 100):
        #    p.camera.view_angle = va
        #    p.write_frame()

        #for el in np.linspace(50, 0, 100):
        #    p.camera.elevation = el
        #    p.write_frame()

        p.close()












