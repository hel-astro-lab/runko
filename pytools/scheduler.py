import sys, os
import pytools  # runko python tools


def str_to_class(classname):
    return getattr(sys.modules[__name__], classname)

class Scheduler:

    def __init__(self, grid):
        self.timer = None
        self.grid = grid

    def is_active_tile(self, tile):
        return True

    def operate(self, op):
    
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
    
            self.grid.send_data(mpid)
            self.grid.recv_data(mpid)
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
    
