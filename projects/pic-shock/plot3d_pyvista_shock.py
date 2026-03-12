import numpy as np
import os
import argparse
import math

import pyvista as pv

import runko
from runko.auto_outdir import resolve_outdir
from runko.mpiio_reader import read_field_snapshot

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



def add_shock_derived(conf):
    """Compute derived charge/mass/B-field quantities on a v5 Configuration.

    Same formulas as pic.py:43-56 so get_normalization() works.
    """
    oppc = 2 * conf.ppc
    omp  = conf.cfl / conf.c_omp
    conf.qe = -(omp**2 * conf.upstream_gamma) / (0.5 * oppc * (1.0 + conf.m0 / conf.m1))
    conf.me = conf.m0 * abs(conf.qe)
    conf.binit = math.sqrt(
        conf.upstream_gamma * oppc * 0.5 * conf.cfl**2
        * conf.me * (1.0 + conf.me / (conf.m1 * abs(conf.qe)))
        * conf.sigma
    )
    conf.stride = 1


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

    return norm




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="3D PyVista shock visualization")
    parser.add_argument("--conf", required=True, help="Path to .ini config file")
    parser.add_argument("-l", "--lap", type=int, required=True, help="Lap number")
    parser.add_argument("-v", "--var", required=True, help="Variable to plot")
    args_cli = parser.parse_args()

    conf = runko.Configuration(args_cli.conf)
    add_shock_derived(conf)
    conf.outdir = resolve_outdir(conf)

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

    lap = args_cli.lap

    #--------------------------------------------------
    # read binary snapshot; fields are (nz, ny, nx) C-order
    fields = read_field_snapshot(os.path.join(conf.outdir, f"flds_{lap}.bin"))

    rho = fields["n0"] + fields["n1"]
    jx  = fields["jx"]
    jy  = fields["jy"]
    jz  = fields["jz"]
    bx  = fields["bx"]
    by  = fields["by"]
    bz  = fields["bz"]

    # normalize
    rho /= get_normalization('rho', conf)

    jx /= get_normalization('jx', conf)
    jy /= get_normalization('jy', conf)
    jz /= get_normalization('jz', conf)

    bx /= get_normalization('bx', conf)
    by /= get_normalization('by', conf)
    bz /= get_normalization('bz', conf)

    dx = conf.stride/conf.c_omp # skindepth resolution
    origin = 0,0,0

    nz, ny, nx = np.shape(rho)
    print(f"field shape (nz,ny,nx) = ({nz},{ny},{nx})")

    #--------------------------------------------------
    # zoom along x (last axis of (nz, ny, nx))
    if False:
        xmin = 0
        xmax = -1
        rho = rho[:, :, xmin:xmax]
        jx  = jx[ :, :, xmin:xmax]
        jy  = jy[ :, :, xmin:xmax]
        jz  = jz[ :, :, xmin:xmax]
        bx  = bx[ :, :, xmin:xmax]
        by  = by[ :, :, xmin:xmax]
        bz  = bz[ :, :, xmin:xmax]

        nz, ny, nx = np.shape(rho)

    #--------------------------------------------------
    # build single PyVista mesh and populate all base fields
    # ravel() on (nz,ny,nx) C-order gives x-fastest, matching VTK's (nx,ny,nz) convention
    Lx = nx*dx
    Ly = ny*dx
    Lz = nz*dx

    midx = 0.5*Lx
    midy = 0.5*Ly
    midz = 0.5*Lz

    Lx_box = 1.0*Lx
    Ly_box = 1.0*Ly

    mesh = pv.ImageData(
        dimensions=(nx, ny, nz),
        spacing=(dx, dx, dx),
        origin=origin,
    )

    mesh['rho']  = rho.ravel()
    mesh['jx']   = jx.ravel()
    mesh['jy']   = jy.ravel()
    mesh['jz']   = jz.ravel()
    mesh['bx']   = bx.ravel()
    mesh['by']   = by.ravel()
    mesh['bz']   = bz.ravel()
    mesh['bvec'] = np.column_stack([bx.ravel(), by.ravel(), bz.ravel()])
    mesh['b2']   = np.sqrt(bx**2 + by**2 + bz**2).ravel()

    #--------------------------------------------------
    p = pv.Plotter(
            #lighting='three lights'
            #off_screen=True,
            )


    #--------------------------------------------------
    #if True: # density
    if var == 'dens': # volume rendering
        clim = (0.0, 4.0)
        print('min/max', np.min(rho), np.max(rho))

        #ind = np.where(rho < 1.5)
        #rho[ind] = 0.0

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
    #if var == 'dens': # streamlines
    if False: # streamlines
        #p.add_mesh(mesh.outline(), color="k")

        b2 = mesh['b2']
        print('min/max', np.min(b2), np.max(b2), np.mean(b2))
        mesh['scalar'] = 1.0 + b2

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

        stream = mesh.streamlines_from_source(
            seed,
            vectors='bvec',
            max_time=1e5,
            #initial_step_length=0.1,
            integration_direction='both',
            #integration_direction='forward',
            )


        #stream, src = mesh.streamlines(
        #    'bvec',
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

        clim = (-0.3, 0.3)

        dbz = bz - 1.0 # -B_0
        dbz[:,:, :520] = 0.0

        print('min/max dbz', np.min(dbz), np.max(dbz))
        mesh['dbz'] = dbz.ravel()

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

        #--------------------------------------------------
        # j scalar field
        #mesh['j'] = (np.sign(jz)*np.sqrt( jx**2 + jy**2 + jz**2 )).ravel()
        mesh['j'] = jy.ravel()

        clim = (-2.0, 2.0)
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
    if False: # current slice

        #--------------------------------------------------
        # j scalar field
        #mesh['j'] = (np.sign(jz)*np.sqrt( jx**2 + jy**2 + jz**2 )).ravel()
        #smoothed = mesh.gaussian_smooth(std_dev=0.7)

        cmap = 'RdBu'
        clim = (-1.5, 1.5)

        # y slice
        if True:
            single_slice = mesh.slice(
                    normal=[0,1,0],
                    origin=(midx,midy,midz),
                    )

            p.add_mesh(single_slice,
                       scalars='jz',
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
                       scalars='jz',
                       cmap=cmap,
                       clim=clim,
                       show_scalar_bar=False,
                       opacity=1.0,
                       )


    #--------------------------------------------------
    # zoom in b field
    if False: # streamlines

        print(mesh)

        #p.add_mesh(mesh.outline(), color="k")

        b2 = mesh['b2']
        print('min/max', np.min(b2), np.max(b2), np.mean(b2))
        mesh['scalar'] = 1.0 + b2

        seed = pv.Plane(center=(midx,midy,midz),
                        i_size=Lx_box,
                        j_size=Ly_box,
                        i_resolution=30,
                        j_resolution=6,
                        direction=(0,0,1),
                        )

        #p.add_mesh(seed.outline(), color='r')

        #stream = mesh.streamlines_from_source(
        #    seed,
        #    vectors='bvec',
        #    max_time=2e4,
        #    initial_step_length=0.01,
        #    integration_direction='both',
        #    )

        stream = mesh.streamlines(
             'bvec',
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


    #--------------------------------------------------
    if False: # rho contours rendering
        print('min/max', np.min(rho), np.max(rho), ' mean', np.mean(rho))

        smoothed = mesh.gaussian_smooth(std_dev=0.3)

        p.add_mesh( smoothed.outline_corners(), color="k")

        contours = smoothed.contour(np.linspace(1.5, 4.0, 10), scalars='rho')
        #contours = contours.smooth_taubin(50)
        p.add_mesh(contours,
                    opacity=0.05,
                    clim=(0.0, 10.0),
                    show_scalar_bar=False,
                    )

        hc  = smoothed.contour([4.0], scalars='rho')
        ##hc = hc.smooth_taubin(10, pass_band=10.0)
        p.add_mesh(hc,
                   opacity=1.0,
                   clim=(0.,10.0),
                   show_scalar_bar=False,
                   )


    #--------------------------------------------------
    if False: # try isocountours rendering
        clim = (0, 10.0)

        rho_masked = rho.copy()
        rho_masked[rho_masked < 2.0] = 0.0

        mesh['rho_masked'] = rho_masked.ravel()
        p.add_mesh( mesh.outline_corners(), color="k" )


        ops = np.zeros(256)
        #ops[0:192] = 0.0
        #ops[192:]  = np.linspace(0.0, 128, 64)

        ops[64:]  = np.linspace(0.0, 256, 192)

        p.add_volume(mesh,
                     scalars='rho_masked',
                     opacity='linear',
                     clim=clim,
                     show_scalar_bar=False,
                     cmap='viridis',
                     shade=False,
                    )



    #--------------------------------------------------
    if False: # try isocountours rendering
        # j scalar field
        jsigned = np.sign(jz)*np.sqrt( jx**2 + jy**2 + jz**2 )
        mesh['jsigned'] = jsigned.ravel()
        #smoothed = mesh.gaussian_smooth(std_dev=0.8)

        p.add_mesh( mesh.outline_corners(), color="k")

        vv = 0.8
        clim = (-1.0, 1.0)

        contours = mesh.contour([0.5*vv, 0.7*vv, 0.9*vv], scalars='jsigned')
        p.add_mesh(contours,
                    opacity=0.3,
                    clim=clim,
                    show_scalar_bar=False,
                    cmap='RdBu',
                    )

        contours = mesh.contour([-0.5*vv, -0.7*vv, -0.9*vv], scalars='jsigned')
        p.add_mesh(contours,
                    opacity=0.3,
                    clim=clim,
                    show_scalar_bar=False,
                    cmap='RdBu',
                    )

        hcp = mesh.contour([+vv], scalars='jsigned')
        p.add_mesh(hcp,
                   opacity=1.0,
                    clim=clim,
                    show_scalar_bar=False,
                    cmap='RdBu',
                   )

        hcm = mesh.contour([-vv], scalars='jsigned')
        p.add_mesh(hcm,
                   opacity=1.0,
                    clim=clim,
                   show_scalar_bar=False,
                   cmap='RdBu',
                   )


    #--------------------------------------------------
    if False: # try isocountours rendering
        # j scalar field
        jabs = np.sqrt(jx**2 + jy**2 + jz**2)
        mesh['jabs'] = jabs.ravel()

        #smoothed = mesh.gaussian_smooth(std_dev=0.8)

        p.add_mesh( mesh.outline_corners(), color="k")

        contours = mesh.contour(np.linspace(0.2, 1.0, 5), scalars='jabs')
        #contours = contours.smooth_taubin(50)
        p.add_mesh(contours,
                    opacity=0.1,
                    clim=(0.3,2.0),
                    show_scalar_bar=False,
                    )

        hc  = mesh.contour([1.0], scalars='jabs')
        #hc = hc.smooth_taubin(10, pass_band=10.0)

        p.add_mesh(hc,
                   opacity=1.0,
                   clim=(0.3,2.0),
                   show_scalar_bar=False,
                   )


        #contours = mesh.contour([1.0], scalars='jabs', method='marching_cubes')
        ##co = grid.contour([0], values, method='flying_edges')
        #dist = np.linalg.norm(mesh.points, axis=1)

        #mesh.plot(
        #        scalars=dist,
        #        smooth_shading=True,
        #        specular=1,
        #        cmap="plasma",
        #        show_scalar_bar=False)

    #--------------------------------------------------
    # center camera on grid; look from top-left-front
    diag = np.sqrt(Lx**2 + Ly**2 + Lz**2)
    cam_dist = 2.0*diag
    cam_dir = np.array([-1.0, -1.0, 1.0])
    cam_dir /= np.linalg.norm(cam_dir)
    cam_pos = (midx + cam_dist*cam_dir[0],
               midy + cam_dist*cam_dir[1],
               midz + cam_dist*cam_dir[2])
    focal   = (midx, midy, midz)
    view_up = (0, 0, 1)
    p.camera_position = [cam_pos, focal, view_up]

    #p.enable_depth_peeling(100)
    #p.enable_anti_aliasing('msaa')
    #p.enable_ssao(kernel_size=32)

    #old
    #p.save_graphic("3d_jz_b_v2.pdf", title="", raster=True, painter=True)

    cpos = p.show(return_cpos=True, auto_close=False)
    #p.camera_position = cpos


    #--------------------------------------------------
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
        p.close()

