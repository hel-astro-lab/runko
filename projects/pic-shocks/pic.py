# -*- coding: utf-8 -*-

from mpi4py import MPI
import numpy as np
import sys, os

# runko + auxiliary modules
import pytools  # runko python tools

# problem specific modules
from init_problem import Configuration_Shocks as Configuration
from init_problem import velocity_profile
from init_problem import density_profile
from init_problem import weigth_profile

# import specific tools from shock toolbox
from injector import MovingInjector
from box_chunking import Chunker
from shock_toolset import ShockToolset


#--------------------------------------------------
rnd_seed_default = 1
np.random.seed(rnd_seed_default)  # global simulation seed

#--------------------------------------------------
# Field initialization
def insert_em_fields(grid, conf, do_initialization=True):

    # into radians
    #btheta = conf.btheta / 180.0 * np.pi
    #bphi = conf.bphi / 180.0 * np.pi
    beta = conf.beta

    for tile in pytools.tiles_all(grid):
        g = tile.get_grids(0)

        if conf.oneD:
            ii,jj,kk = (tile.index, 0, 0)
        elif conf.twoD:
            ii,jj,kk = (*tile.index, 0)
        else:
            ii,jj,kk = tile.index
        
        # insert values into Yee lattices; includes halos from -3 to n+3
        if do_initialization:
            for n in range(-3, conf.NzMesh + 3):
                for m in range(-3, conf.NyMesh + 3):
                    for l in range(-3, conf.NxMesh + 3):
                        if not(conf.use_maxwell_split):

                            # global coordinates
                            #iglob, jglob, kglob = pytools.ind2loc((ii, jj, kk), (l, m, n), conf)
                            #r = np.sqrt(iglob ** 2 + jglob ** 2 + kglob ** 2)

                            # x is parallel direction (along the shock motion)
                            # y is perpendicular to shock front (in 2D)
                            # z is perpendicular to shock front (in 3D)

                            g.bx[l, m, n] = conf.binit*conf.bpar  #np.cos(bphi)
                            g.by[l, m, n] = conf.binit*conf.bplan #conf.binit * bperp #np.sin(bphi) * np.sin(btheta)
                            g.bz[l, m, n] = conf.binit*conf.bperp #np.sin(bphi) * np.cos(btheta)

                            g.ex[l, m, n] = 0.0
                            g.ey[l, m, n] = -beta * g.bz[l, m, n]
                            g.ez[l, m, n] = +beta * g.by[l, m, n]


        # insert damping reference point
        for n in range(conf.NzMesh):
            for m in range(conf.NyMesh):
                for l in range(conf.NxMesh):
                    if not(conf.use_maxwell_split):

                        # copy same setup to reference values that we relax into by calling damp_fields
                        tile.bx_ref[l, m, n] = conf.binit*conf.bpar  #np.cos(bphi)
                        tile.by_ref[l, m, n] = conf.binit*conf.bplan #conf.binit * np.sin(bphi) * np.sin(btheta)
                        tile.bz_ref[l, m, n] = conf.binit*conf.bperp #np.sin(bphi) * np.cos(btheta)

                        tile.ex_ref[l, m, n] = 0.0
                        tile.ey_ref[l, m, n] = -beta * tile.bz_ref[l, m, n]
                        tile.ez_ref[l, m, n] = +beta * tile.by_ref[l, m, n]

    return




# loading of RX damping tiles to the grid.
# NOTE: We need all tiles to have a damping capability because injector is moving 
#       through the grid and will update the damping zone location between `fld1` and `fld2`.
def load_damping_tiles(n, conf):

    for k in range(n.get_Nz()):
        for j in range(n.get_Ny()):
            for i in range(n.get_Nx()):
                # print("{} ({},{}) {} ?= {}".format(n.rank, i,j, n.get_mpi_grid(i,j), ref[j,i]))
		
		
                is_owned = False
                if conf.oneD and n.get_mpi_grid(i) == n.rank():
                    is_owned = True
                    ind = (i,)
                elif conf.twoD and n.get_mpi_grid(i,j) == n.rank():
                    is_owned = True
                    ind = (i,j)
                elif conf.threeD and n.get_mpi_grid(i,j,k) == n.rank():
                    is_owned = True
                    ind = (i,j,k)
                if is_owned:
                    tile = pypic.Tile_wall_RX(conf.NxMesh, conf.NyMesh, conf.NzMesh)

         
                    # add tile
                    pytools.pic.initialize_tile(tile, (i,j,k), n, conf)

                    tile.fld1 = 5
                    tile.fld2 = 6

                    tile.ksupp = 10.0

                    n.add_tile(tile, ind)

    return



from pytools import Scheduler
class ShockScheduler(Scheduler):

    # NOTE we inherit most of the original Scheduler functionality


    # Check if tile is actively used in computations;
    # shock criteria is that it is behind the right wall (injector) and infront 
    # of the chunker. Both locations have a buffer zone of a few tiles to make the
    # check a bit more relaxed
    def is_active_tile(self, tile):
    
        if not(conf.use_injector):
            return True
    
        rightx = tile.mins[0] + self.NxMesh  #right edge of the tile
        leftx  = tile.mins[0]                #left edge of the tile
    
        rwallx = self.rwall.wloc1 +3*self.NxMesh  #right edge of the injector (pad with 3 safety buffer)
        lwallx = self.lwall.walloc -3*self.NxMesh  #left edge of the lwall
    
        if leftx < rwallx and rightx > lwallx:
            return True
        else:
            return False


#-------------------------------------------------- 
#-------------------------------------------------- 
#-------------------------------------------------- 
if __name__ == "__main__":

    # --------------------------------------------------
    # initialize auxiliary tools
    sch  = ShockScheduler() # Set MPI aware task scheduler
    args = pytools.parse_args() # parse command line arguments
    tplt = pytools.TerminalPlot(38,1) # terminal plotting tool

    # create conf object with simulation parameters based on them
    conf = Configuration(args.conf_filename, do_print=sch.is_master)
    sch.conf = conf # remember to update conf to scheduler

    if conf.mpi_task_mode: sch.switch_to_task_mode() # switch scheduler to task mode

    # --------------------------------------------------
    # Timer for profiling
    timer = pytools.Timer()
    timer.start("total")
    timer.start("init")
    timer.do_print = sch.is_example_worker
    timer.verbose = 0  # 0 normal; 1 - debug mode

    sch.timer = timer # remember to update scheduler

    # --------------------------------------------------
    # create output folders
    if sch.is_master: pytools.create_output_folders(conf)

    # --------------------------------------------------
    # load runko
    if conf.threeD:
        # 3D modules
        import pycorgi.threeD as pycorgi   # corgi ++ bindings
        import pyrunko.pic.threeD as pypic # pic c++ bindings
        import pyrunko.emf.threeD as pyfld # fld c++ bindings
    elif conf.twoD:
        # 2D modules
        import pycorgi.twoD as pycorgi     # corgi ++ bindings
        import pyrunko.pic.twoD as pypic   # runko pic c++ bindings
        import pyrunko.emf.twoD as pyfld   # runko fld c++ bindings
    elif conf.oneD:
        # 2D modules
        import pycorgi.oneD as pycorgi     # corgi ++ bindings
        import pyrunko.pic.oneD as pypic   # runko pic c++ bindings
        import pyrunko.emf.oneD as pyfld   # runko fld c++ bindings
    # --------------------------------------------------
    # setup grid
    grid = pycorgi.Grid(conf.Nx, conf.Ny, conf.Nz)
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)
    sch.grid = grid # remember to update scheduler


    # compute initial mpi ranks using Hilbert's curve partitioning
    if conf.use_injector:
        pytools.load_catepillar_track_mpi(grid, conf.mpi_track, conf) 
    else:
        #pytools.load_mpi_y_strides(grid, conf) #equal x stripes

        # advanced Hilbert curve optimization (leave 4 first ranks empty)
        #if grid.size() > 24: # big run
        #    pytools.balance_mpi_3D_rootmem(grid, 4) 
        #else:

        pytools.balance_mpi(grid, conf) #Hilbert curve optimization 


    # load pic tiles into grid
    #pytools.pic.load_tiles(grid, conf)
    load_damping_tiles(grid, conf) # NOTE different from regular intialization

    # --------------------------------------------------
    # simulation restart

    # get current restart file status
    io_stat = pytools.check_for_restart(conf)


    # no restart file; initialize simulation
    if io_stat["do_initialization"]:
        if sch.is_master: print("initializing simulation...")
        lap = 0

        rseed = MPI.COMM_WORLD.Get_rank()
        np.random.seed(rseed + rnd_seed_default)  # sync rnd generator seed for different mpi ranks

        # injecting plasma particles
        if not(conf.use_injector):
            prtcl_stat = pytools.pic.inject(grid, velocity_profile, density_profile, conf)
            if sch.is_example_worker: 
                print("injected:")
                print("     e- prtcls: {}".format(prtcl_stat[0]))
                print("     e+ prtcls: {}".format(prtcl_stat[1]))

        # inserting em grid
        insert_em_fields(grid, conf, do_initialization=True)

        # save a snapshot of the state to disk
        pytools.save_mpi_grid_to_disk(conf.outdir, 0, grid, conf)

    else:
        if sch.is_master:
            print("restarting simulation from lap {}...".format(io_stat["lap"]))

        # insert damping reference values
        insert_em_fields(grid, conf, do_initialization=False)

        # read restart files
        pyfld.read_grids(grid, io_stat["read_lap"], io_stat["read_dir"])
        pypic.read_particles(grid, io_stat["read_lap"], io_stat["read_dir"])

        # set particle types
        for tile in pytools.tiles_all(grid):
            for ispcs in range(conf.Nspecies):
                container = tile.get_container(ispcs)
                container.type = conf.prtcl_types[ispcs] # name container

        # step one step ahead
        lap = io_stat["lap"] + 1

    # --------------------------------------------------
    # static load balancing setup; communicate neighborhood info once

    if sch.is_master: print("load balancing grid..."); sys.stdout.flush()

    # update boundaries
    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    MPI.COMM_WORLD.barrier()

    if sch.is_master: print("loading virtual tiles..."); sys.stdout.flush()

    # load virtual mpi halo tiles
    pytools.pic.load_virtual_tiles(grid, conf)


    # --------------------------------------------------
    # load physics solvers

    if sch.is_master: print("loading solvers..."); sys.stdout.flush()

    sch.fldpropE = pyfld.FDTD2()
    sch.fldpropB = pyfld.FDTD2()

    # enhance numerical speed of light slightly to suppress numerical Cherenkov instability
    sch.fldpropE.corr = conf.c_corr
    sch.fldpropB.corr = conf.c_corr


    # --------------------------------------------------
    # particle puhser

    #sch.pusher = pypic.BorisPusher()
    #sch.pusher = pypic.VayPusher()
    sch.pusher = pypic.HigueraCaryPusher()
    #sch.pusher  = pypic.rGCAPusher()

    #if conf.gammarad > 0:
    #    sch.pusher   = pypic.BorisDragPusher()
    #    sch.pusher.drag = conf.drag_amplitude
    #    sch.pusher.temp = 0.0

    # background field from external pusher components
    if conf.use_maxwell_split:
        pusher.bx_ext = conf.binit*conf.bpar  #np.cos(bphi)
        pusher.by_ext = conf.binit*conf.bplan #conf.binit * bperp #np.sin(bphi) * np.sin(btheta)
        pusher.bz_ext = conf.binit*conf.bperp #np.sin(bphi) * np.cos(btheta)

        pusher.ex_ext = 0.0
        pusher.ey_ext = -conf.beta * conf.binit*conf.bperp
        pusher.ez_ext = +conf.beta * conf.binit*conf.bplan

    # --------------------------------------------------
    # particle interpolator
    sch.fintp = pypic.LinearInterpolator()
    #sch.fintp = pypic.QuadraticInterpolator() #2nd order quadratic
    #sch.fintp = pypic.CubicInterpolator()     #3rd order cubic 3d
    #sch.fintp = pypic.QuarticInterpolator()   #4th order quartic; 3d

    # --------------------------------------------------
    # current deposit
    sch.currint = pypic.ZigZag()
    #sch.currint = pypic.ZigZag_2nd()
    #sch.currint = pypic.ZigZag_3rd()
    #sch.currint = pypic.ZigZag_4th()
    #sch.currint = pypic.Esikerpov_2nd() # 3d only
    #sch.currint = pypic.Esikerpov_4th() # 3d only

    # --------------------------------------------------
    #filter
    sch.flt = pyfld.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)


    # --------------------------------------------------
    # I/O objects
    if sch.is_master: print("loading IO objects..."); sys.stdout.flush()

    # quick field snapshots
    #fld_writer = pyfld.MasterFieldsWriter(
    fld_writer = pyfld.FieldsWriter(
        conf.outdir,
        conf.Nx,
        conf.NxMesh,
        conf.Ny,
        conf.NyMesh,
        conf.Nz,
        conf.NzMesh,
        conf.stride,
    )


    # tracked particles; only works with no injector (id's are messed up when coming out from injector)
    prtcl_writers = []

    #for ispc in [0, 1, 2]: #electrons & positrons
    for ispc in [0]: #electrons
        mpi_comm_size = grid.size() if not(conf.mpi_task_mode) else grid.size() - 1
        n_local_tiles = int(conf.Nx*conf.Ny*conf.Nz/mpi_comm_size)
        #print("average n_local_tiles:", n_local_tiles, " / ", mpi_comm_size)

        prtcl_writer = pypic.TestPrtclWriter(
                conf.outdir,
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh,
                conf.ppc,
                n_local_tiles, #len(grid.get_local_tiles()),
                conf.n_test_prtcls,)
        prtcl_writer.ispc = ispc
        prtcl_writers.append(prtcl_writer)


    # momens of particle distribution
    #mom_writer = pypic.MasterPicMomentsWriter(
    mom_writer = pypic.PicMomentsWriter(
        conf.outdir,
        conf.Nx,
        conf.NxMesh,
        conf.Ny,
        conf.NyMesh,
        conf.Nz,
        conf.NzMesh,
        conf.stride_mom,
    )

    # 3D box peripherals
    if conf.threeD:
        st = 1 # stride
        slice_xy_writer = pyfld.FieldSliceWriter( conf.outdir, 
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh, st, 0, 1)
        slice_xz_writer = pyfld.FieldSliceWriter( conf.outdir, 
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh, st, 1, 1)
        slice_yz_writer = pyfld.FieldSliceWriter( conf.outdir, 
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh, st, 2, 1)

        # location of the slice (grid index)
        #slice_xy_writer.ind = int(0.0*conf.Lz) # bottom slice
        #slice_xz_writer.ind = int(0.0*conf.Ly) # side wall 1
        #slice_yz_writer.ind = int(0.0*conf.Lx) # side wall 2



    # --------------------------------------------------
    # shock specific tools
    if sch.is_master: print("loading shock utilities..."); sys.stdout.flush()

    sch.NxMesh = conf.NxMesh # needed for calculating the buffer zone 

    shock_toolset = ShockToolset(conf) # shock specific IO tools (i.e., tracking of shock front)

    sch.lwall = pypic.Piston() # reflecting leftmost wall

    if conf.use_maxwell_split: # remember to reset reflectors EM fields if split is used
        sch.lwall.bxwall = 0.0
        sch.lwall.bywall = 0.0
        sch.lwall.bzwall = 0.0 # use background magnetic field that is set by pusher

        sch.lwall.exwall = 0.0
        sch.lwall.eywall = +conf.beta*conf.binit*conf.bperp # cancel ey of pusher so ey = 0
        sch.lwall.ezwall = -conf.beta*conf.binit*conf.bplan # cancel ez of pusher so ez = 0


    # set piston wall speed (for standard reflector it is non-moving so gam = 0)
    sch.lwall.walloc = 15.0  # leave min of >10 cell spacing between the wall for boundary conditions

    if conf.wallgamma > 0.0: #moving wall
        sch.lwall.gammawall = conf.wallgamma
        sch.lwall.betawall = np.sqrt(1.0 - 1.0 / conf.wallgamma ** 2.0)
    else: #stationary wall
        sch.lwall.betawall = 0.0
        sch.lwall.gammawall = 1.0

    # --------------------------------------------------
    # moving injector
    sch.rwall = MovingInjector(conf)

    # inject initial stripe of prtcls
    if io_stat["do_initialization"]:

        if conf.use_injector:
            if sch.is_master:
                print("injecting prtcls with moving injector")
                sys.stdout.flush()

            prtcl_stat = sch.rwall.inject(grid, lap, velocity_profile, density_profile, conf)
            if sch.is_master:
                print('injected e^-:', prtcl_stat[0])
                print('injected e^+:', prtcl_stat[1])
        sch.rwall.step(lap)


    # --------------------------------------------------
    # box slider
    slider = Chunker(conf)

    # restaring; update internal counters to match lap
    if not(io_stat["do_initialization"]):

        # NOTE: we perform this by bruteforce re-looping over steps
        #if conf.use_injector:
        sch.rwall.step(0)

        for plap in range(0, lap):
            do_slide = slider.slide(plap, grid, conf)

            if do_slide: sch.lwall.walloc = slider.xmin  # set reflecting box to end of the box

            sch.rwall.step(plap)

            shock_toolset.wloc0 = slider.xmin #downstream chunking sets left box limit
            shock_toolset.wloc1 = sch.rwall.wloc1 #upstream injector sets right box limit

        if sch.is_master:
            print("lwall walloc", sch.lwall.walloc)
            print("shocktl shock_loc", shock_toolset.shock_loc)
            print("rwall wloc0", sch.rwall.wloc0)
            print("rwall wloc1", sch.rwall.wloc1)
            print("slider xmin", slider.xmin)

    # --------------------------------------------------
    # --------------------------------------------------
    # --------------------------------------------------
    # end of initialization
    sch.timer.stop("init")
    sch.timer.stats("init")

    if sch.is_master: print('init: starting simulation...'); sys.stdout.flush()


    ##################################################
    # simulation time step loop

    # simulation loop
    time = lap * (conf.cfl / conf.c_omp)
    for lap in range(lap, conf.Nt + 1):

        # --------------------------------------------------
        # sliding box 
        t1 = sch.timer.start_comp("sliding_box")
        do_slide = slider.slide(lap, grid, conf)
        if do_slide:
            sch.lwall.walloc = slider.xmin  # set reflecting box to the end of the box
        sch.timer.stop_comp(t1)

        sch.rwall.step(lap)

        # moving right injector
        t1 = sch.timer.start_comp("rwall_bc")
        sch.rwall.damp_em_fields(grid, lap, conf)
        sch.timer.stop_comp(t1)

        if conf.use_injector:
            t1 = sch.timer.start_comp("rwall")
            prtcl_stat = sch.rwall.inject(grid, lap, velocity_profile, density_profile, conf)
            sch.timer.stop_comp(t1)

        # --------------------------------------------------
        # comm E and B
        sch.operate( dict(name='mpi_b0', solver='mpi', method='b', ) )
        sch.operate( dict(name='mpi_e0', solver='mpi', method='e', ) )
        sch.operate( dict(name='upd_bc', solver='tile',method='update_boundaries', args=[grid,[1,2] ], nhood='local', ) )

        # --------------------------------------------------
        # push B half
        sch.operate( dict(name='push_half_b1', solver='fldpropB', method='push_half_b', nhood='local',) )
        sch.operate( dict(name='wall_bc',      solver='lwall',    method='field_bc',    nhood='local',) )

        # comm B
        sch.operate( dict(name='mpi_b1',  solver='mpi', method='b',                 ) )
        sch.operate( dict(name='upd_bc ', solver='tile',method='update_boundaries',args=[grid, [2,] ], nhood='local',) )

        # --------------------------------------------------
        # move particles (only locals tiles)

        # interpolate fields and push particles in x and u
        sch.operate( dict(name='interp_em', solver='fintp',  method='solve', nhood='local', ) )
        sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', args=[0]) ) # e^-
        sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', args=[1]) ) # e^+


        # clear currents; need to call this before wall operations since they can deposit currents too 
        sch.operate( dict(name='clear_cur', solver='tile',   method='clear_current', nhood='all', ) )

        # apply moving/reflecting walls
        sch.operate( dict(name='walls',     solver='lwall', method='solve', nhood='local', ) )


        # --------------------------------------------------
        # advance half B 
        sch.operate( dict(name='push_half_b2', solver='fldpropB', method='push_half_b', nhood='local', ) )
        sch.operate( dict(name='wall_bc',      solver='lwall',    method='field_bc',    nhood='local', ) )

        # comm B
        sch.operate( dict(name='mpi_b2', solver='mpi', method='b',                 ) )
        sch.operate( dict(name='upd_bc', solver='tile',method='update_boundaries', args=[grid, [2,] ], nhood='local',) )


        # --------------------------------------------------
        # push E
        sch.operate( dict(name='push_e',   solver='fldpropE', method='push_e',  nhood='local', ) )
        sch.operate( dict(name='wall_bc',  solver='lwall',    method='field_bc',nhood='local', ) )

        # --------------------------------------------------
        # particle communication (only local/boundary tiles)

        # local and global particle exchange 
        sch.operate( dict(name='check_outg_prtcls',     solver='tile',  method='check_outgoing_particles',     nhood='local', ) )
        sch.operate( dict(name='pack_outg_prtcls',      solver='tile',  method='pack_outgoing_particles',      nhood='boundary', ) )

        sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p1',                           nhood='all', ) )
        sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p2',                           nhood='all', ) )

        sch.operate( dict(name='unpack_vir_prtcls',     solver='tile',  method='unpack_incoming_particles',    nhood='virtual', ) )
        sch.operate( dict(name='check_outg_vir_prtcls', solver='tile',  method='check_outgoing_particles',     nhood='virtual', ) )
        sch.operate( dict(name='get_inc_prtcls',        solver='tile',  method='get_incoming_particles',       nhood='local', args=[grid,]) )

        sch.operate( dict(name='del_trnsfrd_prtcls',    solver='tile',  method='delete_transferred_particles', nhood='local', ) )
        sch.operate( dict(name='del_vir_prtcls',        solver='tile',  method='delete_all_particles',         nhood='virtual', ) )

        # --------------------------------------------------
        # current calculation; charge conserving current deposition
        # clear virtual current arrays for boundary addition after mpi, send currents, and exchange between tiles
        sch.operate( dict(name='comp_curr',     solver='currint', method='solve',             nhood='local', ) )
        sch.operate( dict(name='clear_vir_cur', solver='tile',    method='clear_current',     nhood='virtual', ) )
        sch.operate( dict(name='mpi_cur',       solver='mpi',     method='j',                 nhood='all', ) )
        sch.operate( dict(name='cur_exchange',  solver='tile',    method='exchange_currents', nhood='local', args=[grid,], ) )

        # --------------------------------------------------
        # filter
        for fj in range(conf.npasses):

            # flt uses halo=2 padding so only every 3rd (0,1,2) pass needs update
            if fj % 2 == 0:
                sch.operate( dict(name='mpi_cur_flt', solver='mpi', method='j', ) )
                sch.operate( dict(name='upd_bc',      solver='tile',method='update_boundaries',args=[grid, [0,] ], nhood='local', ) )

                #TODO is wall_bc needed here?

                MPI.COMM_WORLD.barrier()
            sch.operate( dict(name='filter', solver='flt', method='solve', nhood='local', ) )


        # --------------------------------------------------
        # add antenna contribution for exciting waves
        #sch.antenna.update_rnd_phases()
        #antenna.get_brms(grid) # debug tracking
        #sch.operate( dict(name='add_antenna', solver='antenna', method='add_ext_cur', nhood='local', ) )

        # --------------------------------------------------
        # add current to E
        sch.operate( dict(name='add_cur', solver='tile', method='deposit_current', nhood='local', ) )
        sch.operate( dict(name='wall_bc', solver='lwall', method='field_bc', nhood='local', ) )


        ##################################################
        # data reduction and I/O

        sch.timer.lap("step")
        if lap % conf.interval == 0:
            if sch.is_master:
                print("--------------------------------------------------")
                print("------ lap: {} / t: {}".format(lap, time))

            sch.timer.stats("step")
            sch.timer.comp_stats()
            sch.timer.purge_comps()

            # io/analyze (independent)
            sch.timer.start("io")
            # shrink particle arrays
            for tile in pytools.tiles_all(grid):
                tile.shrink_to_fit_all_particles()

            # barrier for quick writers
            MPI.COMM_WORLD.barrier()

            # shallow IO
            # NOTE: do moms before other IOs to keep rho field up-to-date
            mom_writer.write(grid, lap)  # pic distribution moments; 
            fld_writer.write(grid, lap)  # quick field snapshots

            for pw in prtcl_writers:
                pw.write(grid, lap)  # particle tracking
            
            #pytools.save_mpi_grid_to_disk(conf.outdir, lap, grid, conf) # MPI grid

            #box peripheries 
            if conf.threeD:
                slice_xy_writer.write(grid, lap)
                slice_xz_writer.write(grid, lap)
                slice_yz_writer.write(grid, lap)


            #--------------------------------------------------
            # shock analysis toolset

            # set box limits
            #lwall.walloc = slider.xmin       # set reflecting box to end of the box
            shock_toolset.wloc0 = slider.xmin #downstream chunking sets left box limit
            shock_toolset.wloc1 = sch.rwall.wloc1 #upstream injector sets right box limit

            shock_toolset.get_density_profile(grid, lap, conf)
            shock_toolset.find_shock_front(grid, lap, conf, xmin=slider.xmin)

            shock_toolset.read_fields(grid, lap, conf)
            shock_toolset.read_prtcls(grid, lap, conf)

            shock_toolset.save(grid, lap, conf)


            #--------------------------------------------------
            # terminal plot 
            if sch.is_master:
                tplt.col_mode = False # use terminal color / use ASCII art

                if conf.twoD:
                    tplt.plot_panels( (2,1),
                    dict(axs=(0,0), data=fld_writer.get_slice(0)/conf.p_norm, name='ne', cmap='viridis', vmin=0,  vmax=2),
                    dict(axs=(1,0), data=fld_writer.get_slice(6)/conf.j_norm, name='jx', cmap='RdBu',    vmin=-1, vmax=1),
                                     )

                elif conf.threeD:
                    tplt.plot_panels( (1,1),
                    dict(axs=(0,0), data=slice_xy_writer.get_slice(8)/conf.j_norm ,   name='jx (xy)', cmap='RdBu'   ,vmin=-2, vmax=2),
                    #dict(axs=(0,0), data=slice_xy_writer.get_slice(9)/conf.p_norm   , name='ne (xy)', cmap='viridis',vmin= 0, vmax=4),
                    #dict(axs=(1,0), data=slice_xy_writer.get_slice(3)/conf.b_norm ,   name='bx (xy)', cmap='RdBu'   ,vmin=-2, vmax=2),
                    #dict(axs=(1,1), data=slice_xy_writer.get_slice(5)/conf.b_norm ,   name='bz (xy)', cmap='RdBu'   ,vmin=-2, vmax=2),
                    )

            #--------------------------------------------------
            #print statistics
            if sch.is_master:
                print('simulation time    {:7d} ({:7.1f} omp)   {:5.1f}%'.format( int(lap), time, 100.0*lap/conf.Nt))
                print('reflecting wall at {:7d} ({:7.1f} c/omp) {:5.1f}%'.format( int(sch.lwall.walloc), sch.lwall.walloc/conf.c_omp, 100.0*sch.lwall.walloc/slider.Lx))
                print('shock front at     {:7d} ({:7.1f} c/omp) {:5.1f}%'.format( int(shock_toolset.shock_loc), shock_toolset.shock_loc/conf.c_omp, 100.0*shock_toolset.shock_loc/slider.Lx))
                print('injector at        {:7d} ({:7.1f} c/omp) {:5.1f}%'.format( int(sch.rwall.wloc1), sch.rwall.wloc1/conf.c_omp, 100.0*sch.rwall.wloc1/slider.Lx))

            #--------------------------------------------------
            # deep IO
            if conf.full_interval > 0 and (lap % conf.full_interval == 0) and (lap > 0):
                pyfld.write_grids(grid, lap, conf.outdir + "/full_output/")
                pypic.write_particles(grid, lap, conf.outdir + "/full_output/")


            # restart IO (overwrites)
            if (lap % conf.restart == 0) and (lap > 0):

                # flip between two sets of files
                io_stat["deep_io_switch"] = 1 if io_stat["deep_io_switch"] == 0 else 0

                pyfld.write_grids(
                    grid, io_stat["deep_io_switch"] + io_stat['restart_num'],
                    conf.outdir + "/restart/"
                )

                pypic.write_particles(
                    grid, io_stat["deep_io_switch"] + io_stat['restart_num'],
                    conf.outdir + "/restart/"
                )

                # if successful adjust info file
                MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write
                if grid.rank() == 0:
                    with open(conf.outdir + "/restart/laps.txt", "a") as lapfile:
                        lapfile.write("{},{}\n".format(
                            lap, 
                            io_stat["deep_io_switch"]+io_stat['restart_num']))

            MPI.COMM_WORLD.barrier() # extra barrier to synch everybody after IOs

            timer.stop("io")
            timer.stats("io")
            #--------------------------------------------------

            timer.start("step")  # refresh lap counter (avoids IO profiling)

            sys.stdout.flush()

        # next step
        time += conf.cfl / conf.c_omp
        #MPI.COMM_WORLD.barrier() # extra barrier to synch everybody
    # end of loop


    # --------------------------------------------------
    # end of simulation

    timer.stop("total")
    timer.stats("total")

