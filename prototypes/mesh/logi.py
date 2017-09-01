import numpy as np
import math
import copy

import conf 


##################################################
# Auxiliary functions
def xwrap(i):
    while i < 0:
        i += conf.Nx
    while i >= conf.Nx:
        i -= conf.Nx
    return i

def ywrap(j):
    while j < 0:
        j += conf.Ny
    while j >= conf.Ny:
        j -= conf.Ny
    return j


#hard boundaries
#def xwrap(i):
#    while i < 0:
#        return None
#    while i >= conf.Nx:
#        return None
#    return i
#
#def ywrap(j):
#    while j < 0:
#        return None
#    while j >= conf.Ny:
#        return None
#    return j


# Copy & append
def cappend(arr, x):
    tmpa = copy.deepcopy( arr )
    tmpx = copy.deepcopy( x )
    tmpa.append( tmpx )
    return tmpa

# copy delete
def cdel(arr, i):
    tmp = copy.deepcopy( arr )
    del tmp[i]
    return tmp


##################################################

# spatial grid located inside one cell
class Grid:
    xmin = 0.0
    xmin = 1.0
    ymin = 0.0
    ymax = 1.0

    def __init__(self, xmin, xmax, ymin, ymax, NxCell, NyCell):
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        self.mesh = np.zeros( (NxCell, NyCell) )

        #self.mesh[0,0] = 1.0
        #self.mesh[1,0] = 2.0
        #self.mesh[0,1] = 3.0
        #self.mesh[1,1] = 4.0


    



#Main computational cell 
class cell:


    grid = None

    state = 'active'
    cid   = 0
    owner = 0

    communications              = None
    number_of_virtual_neighbors = None
    top_owner                   = None


    i = None
    j = None

    def __init__(self, i, j, owner):
        self.i = i
        self.j = j
        self.owner = owner

    def index(self):
        return ( self.i, self.j )

    #relative neighbors
    def neighs(self, ir, jr):
        i = xwrap(self.i + ir)
        j = ywrap(self.j + jr)
        return (i, j)

    def full_neighborhood(self):
        nb = []
        for ir in [-1, 0, 1]:
            for jr in [-1, 0, 1]:
                if not( (ir == 0) and (jr == 0) ):
                    nb.append( self.neighs(ir, jr) )
        return nb


##################################################
#Main computational node holding cells and 
# dealing with inter-cell communications 
class node:

    cells    = []
    virtuals = []

    send_queue = 'building'
    send_queue_cells   = []
    send_queue_address = []

    adopted_index  = []
    adopted_parent = []
    kidnap_index   = []
    kidnapper      = []

    purged = []

    total_work = None
    work_goal  = None



    def __init__(self, rank, Nx, Ny):
        self.rank = rank
        self.mpiGrid = np.zeros( (Nx, Ny) )
        

    #Purge previous id list and re-create it
    # we also re-tag every cell during the process
    #def numerize(self):
    #    self.cellList = range(len(self.cells))
    #    for cid, c in enumerate(self.cells):
    #        c.owner = self.rank
    #        c.cid   = cid

    #def indexify(self):
    #    self.indexes = []
    #    for c in self.cells:
    #        self.indexes.append( c.index() )

    #Simple ownership test using cell's owner id
    #def is_mine(self, c):
    #    if c.owner == self.rank:
    #        return True
    #    else:
    #        return False

    #Check if given cell id is in our cell list 
    # this is more rigorous test of ownership
    def is_local(self, indx):
        local = False
        for c in self.cells:
            if c.index() == indx:
                local = True
                break
        return local

    
    # Get ALL virtual neighbors relevant for the current node
    def get_all_virtuals(self):
        neighbors = []
        locals    = []
        virtuals  = []

        for c in self.cells:
            locals.append( c.neighs(0, 0) )
            neighbors.extend( c.full_neighborhood() )

        for (i, j) in neighbors:
            if i == None or j == None:
                #boundary value that does not exist, i.e., edge 
                break

            for (il, jl) in locals:
                if (i == il) and (j == jl):
                    #neighbor is local so we omit it
                    break
            else:
                for (iv, jv) in virtuals:
                    if (iv == i) and (jv == j):
                        #neighbor is already on the list so we omit it
                        break
                else:
                    # unique virtual index
                    virtuals.append( (i,j) )

        return virtuals
        
    # returns index of the neighboring cell
    # given in relative coordinates w.r.t. c
    def get_neighbor_index(self, c, i, j):
        return c.neighs(i, j)

    # returns a pointer to the cell given with
    # relative coordinates w.r.t. c
    def get_neighbor_cell(self, c, i, j):
        #get global index
        indx = self.get_neighbor_index(c, i, j)

        if self.is_local(indx):
            #loop over all cells to find it
            for q, c in enumerate(self.cells):
                if c.index() == indx:
                    return c
        else: #virtual
            return -1

    #number of virtual (non-local) neighboring 
    # cells of the given cell
    def number_of_virtual_neighbors(self, c):
        neigs = c.full_neighborhood()
        N = 0
        for indx in neigs:
            if (indx[0] == None or indx[1] == None):
                continue
            if not( self.is_local(indx) ):
                N += 1
        return N

    #owners of the virtual cells around cell c
    def virtual_neighborhood(self, c):
        neigs = c.full_neighborhood()
        owners = []
        for indx in neigs:
            if (indx[0] == None or indx[1] == None):
                continue
            if not( self.is_local(indx) ):
                owners.append( int(self.mpiGrid[indx]) )

        #return unique list of owners, and the mode of the array
        return np.unique(owners).tolist(), np.argmax(np.bincount(owners))



    def clear_virtuals(self):
        self.virtuals = []

    # pack boundary cells to be communicated to neighbors
    def pack_all_virtuals(self):
        self.send_queue         = 'building'
        self.send_queue_cells   = []
        self.send_queue_address = []

        packed = []
        for c in self.cells:
            #if self.rank == 1:
            #    print self.rank, " c:", c.index()
            N = self.number_of_virtual_neighbors(c)
            if N > 0:
                owners, top_owner = self.virtual_neighborhood(c)
                #if self.rank == 1:
                #    print self.rank, " packing", c.index(), " for ", owners

                c.top_owner      = top_owner 
                c.communications = len(owners)
                c.number_of_virtual_neighbors = N

                if not(c.index() in packed):
                    self.send_queue_cells = cappend(self.send_queue_cells, c)
                    self.send_queue_address = cappend(self.send_queue_address, owners)
                    packed = cappend( packed, c.index() )
                #else:
                    #if self.rank == 1:
                    #    print "not packing", c.index()," to ", packed

        self.send_queue = 'ready'

        #if self.rank == 1:
        #    for i, c in enumerate(self.send_queue_cells):
        #        print "in queue: ", c.index(), " for ",self.send_queue_address[i]


    def clear_queue(self):
        self.send_queue         = 'cleared'
        self.send_queue_cells   = []
        self.send_queue_address = []


    def rank_virtuals(self):
        Nvirs      = []
        Ncoms      = []
        Tvirs      = []
        owners     = []
        index_list = []

        self.adopted_index  = []
        self.adopted_parent = []

        for c in self.virtuals:
            Nvirs.append( c.number_of_virtual_neighbors )
            Ncoms.append( c.communications )
            Tvirs.append( c.top_owner )
            owners.append( c.owner )
            index_list.append( c.index() )


        return Nvirs, Ncoms, Tvirs, owners, index_list

        #printing
        #indx_sorted = np.argsort(Nvirs).tolist()

        #for i in indx_sorted:
            #print "virtual cell ({},{}): {} with # neighborhood {} and owner {}".format(
            #        index_list[i][1],
            #        index_list[i][0],
            #        Nvir[i], 
            #        Ncom[i],
            #        owners[i]
            #        )

        #adoption
        #adopt last from virtual to real cells
        #last_indx = indx_sorted[-1]

        #self.adopted_index.append( self.virtuals[last_indx].index() )
        #self.adopted_parent.append( self.virtuals[last_indx].owner )



        #numpy list method
        #self.virtuals[last_indx].owner = self.rank #mark new owner
        #self.cells    = np.append( self.cells, np.copy(self.virtuals[last_indx]) )
        #self.virtuals = np.delete( self.virtuals, last_indx )

        #deepcopy method
        #self.virtuals[last_indx].owner = self.rank #mark new owner
        #self.cells    = cappend( self.cells, self.virtuals[last_indx] )
        #self.virtuals = cdel( self.virtuals, last_indx )

        #python list method
        #vc = self.virtuals[ last_indx ]
        #vc.owner = self.rank
        #self.cells.append( vc )





    def adopt(self):
        for indx in self.adopted_index:
            print "node {} got adopt command for ({},{})".format(self.rank, indx[0],indx[1])

            for c in self.virtuals:
                if c.index() == indx:
                    c_tmp       = copy.deepcopy( c )
                    c_tmp.owner = self.rank

                    (i,j) = indx
                    self.mpiGrid[i,j] = self.rank

                    print "   ...adopting ({},{}) to be {}".format(c_tmp.index()[0],c_tmp.index()[1], c_tmp.owner)
                    self.cells  = cappend( self.cells, c_tmp )
                    break


            #loop over whole virtual list in case we have multiples
            to_be_purged = []
            for i, c in enumerate( self.virtuals ):
                if c.index() == indx:
                    to_be_purged.append( i )

            print "Now I'll remove these virtuals", to_be_purged
            for i in reversed( to_be_purged ):
                    self.virtuals = cdel( self.virtuals, i )



    def purge(self):
        indxs = []
        for indx in self.purged:
            for i, c in enumerate(self.cells):
                if c.index() == indx:
                    indxs.append( i )
                    self.cells = cdel( self.cells, i )
                    break

        self.purged = []
        #self.purged = np.array([])

