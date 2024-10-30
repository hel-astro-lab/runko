import sys, os
import pytools  # runko python tools
from mpi4py import MPI



def str_to_class(classname):
    return getattr(sys.modules[__name__], classname)

class Scheduler:

    def __init__(self, print_banner=True):
        self.timer = None
        self.grid  = None

        self.rank = MPI.COMM_WORLD.Get_rank() 
        self.mpi_comm_size = MPI.COMM_WORLD.Get_size() 

        # default mode where root prints
        self.is_master         = True if self.rank == 0 else False # root rank
        self.is_example_worker = True if self.rank == 0 else False # example work rank

        if self.is_master and print_banner:
            pytools.print_banner()
            print("sch : running runko with {} MPI rank(s)".format(self.mpi_comm_size))

        self.mpi_task_mode = False
        self.master_rank       = 0
        self.example_work_rank = 0

        self.debug = False # debug mode

    # swithc from all-in mode to task mode
    def switch_to_task_mode(self,):
        self.mpi_task_mode = True

        if self.mpi_comm_size > 1:
            self.example_work_rank = 1

        if self.is_master:
            print("sch : operating in task mode; rank {} is master; {} is example worker".format(
                self.master_rank, self.example_work_rank))

        # set master rank
        if self.rank == self.master_rank: 
            self.is_master = True
        else:
            self.is_master = False

        # set example worker
        if self.rank == self.example_work_rank: 
            self.is_example_worker = True
        else:
            self.is_example_worker = False


    def is_active_tile(self, tile):
        return True

    def operate(self, op):

        if self.debug: # additional debug printing
            print('R:', self.rank, op) # debug print
            MPI.COMM_WORLD.Barrier() 
            if self.is_master:
                print('',sys.stdout.flush())
    
        #--------------------------------------------------
        # default values
        if not('nhood' in op):
            op['nhood'] = 'all'
    
        if not('args' in op):
            op['args'] = []
    
        #--------------------------------------------------
        # parse neighborhood
        non_boundary = False
    
        if op['nhood'] == 'all':
            tile_iterator = pytools.tiles_all
            non_boundary = True
    
        elif op['nhood'] == 'local':
            non_boundary = True
            tile_iterator = pytools.tiles_local
    
        elif op['nhood'] == 'virtual':
            tile_iterator = pytools.tiles_virtual
        elif op['nhood'] == 'boundary':
            tile_iterator = pytools.tiles_boundary
    
        #-------------------------------------------------- 
        # operate directly to tile 
        if op['solver'] == 'tile':
    
            # actual loop
            t1 = self.timer.start_comp(op['name'])
            for tile in tile_iterator(self.grid):
    
                # skip-non active non-boundary tiles
                if non_boundary:
                    is_active = self.is_active_tile(tile)
                    if not(is_active):
                        continue
    
                method = getattr(tile, op['method'])
                method(*op['args'])
    
            self.timer.stop_comp(t1)
    
        #-------------------------------------------------- 
        # MPI
        elif op['solver'] == 'mpi':
    
            if op['method'] == 'j': mpid = 0
            if op['method'] == 'e': mpid = 1
            if op['method'] == 'b': mpid = 2
            if op['method'] == 'p1': mpid = 3
            if op['method'] == 'p2': mpid = 4
    
            t1 = self.timer.start_comp(op['name'])
    
            self.grid.recv_data(mpid)
            self.grid.send_data(mpid)
            self.grid.wait_data(mpid)
    
            self.timer.stop_comp(t1)
    
        #-------------------------------------------------- 
        #normal solver
        else:
            #solver = str_to_class(op['solver'])
            solver = getattr(self, op['solver'])
            method = getattr(solver, op['method'])
    
            # actual loop
            t1 = self.timer.start_comp(op['name'])
    
            for tile in tile_iterator(self.grid):
    
                # skip-non active non-boundary tiles
                if non_boundary:
                    is_active = self.is_active_tile(tile)
                    if not(is_active):
                        continue
    
                single_args = [tile] + op['args']
                method(*single_args)
    
            self.timer.stop_comp(t1)
    
