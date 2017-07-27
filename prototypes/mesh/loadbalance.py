import numpy as np
import math
from pylab import *
import scipy
import copy

import palettable as pal
from matplotlib import cm

#cmap = cm.get_cmap('inferno_r')
cmap = pal.wesanderson.Moonrise1_5.mpl_colormap





##################################################
#Auxiliary functions to help with the grid

def imshow(ax, grid):

    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)

    extent = [xmin, xmax, ymin, ymax]

    mgrid = np.ma.masked_where(grid == -1.0, grid)
    ax.imshow(mgrid,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = cmap,
              vmin = 0.0,
              vmax = Nrank-1,
              #vmax = Nrank,
              #alpha=0.2
              )


def xwrap(i):
    while i < 0:
        i += Nx
    while i >= Nx:
        i -= Nx
    return i

def ywrap(j):
    while j < 0:
        j += Ny
    while j >= Ny:
        j -= Ny
    return j



def plot_grid(ax, nodes):
    tmp_grid = np.ones( (Nx, Ny) ) * -1.0

    for n in nodes:
        for c in n.cells:
            (i, j) = c.index()
            tmp_grid[i,j] = n.rank
            #tmp_grid[i,j] += 1

    imshow(ax, tmp_grid)




def plot_node(ax, n):
    tmp_grid = np.ones( (Nx, Ny) ) * -1.0

    for c in n.cells:
        (i, j) = c.index()
        #tmp_grid[i,j] = 1.0
        tmp_grid[i,j] = c.owner


    #virtuals = n.get_all_virtuals()
    #for (i,j) in virtuals:
    #    tmp_grid[i,j] += 0.5

    for c in n.virtuals:
        (i,j) = c.index()
    #    #tmp_grid[i,j] += 0.5
        tmp_grid[i,j] = c.owner


    imshow(ax, tmp_grid)

    for c in n.cells:
        (i, j) = c.index()
        ix = xmin + xmax*(i+0.5)/Nx
        jy = ymin + ymax*(j+0.5)/Ny

        Nv = n.number_of_virtual_neighbors(c)
        label = str(Nv)
        ax.text(jy, ix, label, ha='center',va='center')





##################################################

#Main computational cell 
class cell:

    data=0

    state = 'active'
    cid   = 0
    owner = 0

    communications = 0

    i = 0
    j = 0

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



#Main computational node holding cells and 
# dealing with inter-cell communications 
class node:

    #cells    = np.array([], dtype=np.object) #, copy=True)
    #virtuals = np.array([], dtype=np.object) #, copy=True)
    cells    = []
    virtuals = []


    send_queue = 'building'
    send_queue_cells   = []
    send_queue_address = []

    adopted_index = []
    adopted_parent = []

    purged = []
    #purged = np.array([])


    def __init__(self, rank):
        self.rank = rank
        

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
            if not( self.is_local(indx) ):
                N += 1
        return N

    #owners of the virtual cells around cell c
    def virtual_neighborhood(self, c):
        neigs = c.full_neighborhood()
        owners = []
        for indx in neigs:
            if not( self.is_local(indx) ):
                owners.append( np.int(grid[indx]) )
        return np.unique(owners)

    def clear_virtuals(self):
        self.virtuals = []

    # pack boundary cells to be communicated to neighbors
    def pack_virtuals(self):
        self.send_queue         = 'building'
        self.send_queue_cells   = []
        self.send_queue_address = []

        packed = []
        for c in self.cells:
            N = self.number_of_virtual_neighbors(c)
            if N > 0:
                owners = self.virtual_neighborhood(c)

                c.communications = len(owners)
                c.number_of_virtual_neighbors = N

                #self.send_queue_cells.append(c)
                #self.send_queue_address.append(owners)
                if not(c.index() in packed):
                    self.send_queue_cells = cappend(self.send_queue_cells, c)
                    self.send_queue_address = cappend(self.send_queue_address, owners)
                    packed = cappend( packed, c.index() )


        self.send_queue = 'ready'

    def clear_queue(self):
        self.send_queue         = 'cleared'
        self.send_queue_cells   = []
        self.send_queue_address = []


    def rank_virtuals(self):
        Nvir = []
        Ncom = []
        owners = []
        index_list = []

        self.adopted_index = []
        self.adopted_parent = []

        for c in self.virtuals:
            Nvir.append( c.number_of_virtual_neighbors )
            Ncom.append( c.communications )
            owners.append( c.owner )
            index_list.append( c.index() )
        indx_sorted = np.argsort(Nvir)
        
        #for i in indx_sorted:
            #print "virtual cell ({},{}): {} with # neighborhood {} and owner {}".format(
            #        index_list[i][1],
            #        index_list[i][0],
            #        Nvir[i], 
            #        Ncom[i],
            #        owners[i]
            #        )

        #adopt last from virtual to real cells
        last_indx = indx_sorted[-1]

        self.adopted_index.append( self.virtuals[last_indx].index() )
        self.adopted_parent.append( self.virtuals[last_indx].owner )



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





    def adopt(self, rindx):

        print "node {} got adopt command for ({},{})".format(self.rank, rindx[0],rindx[1])

        #print "to be loopped:", self.adopted_index
        for c in self.virtuals:
            #print "...virc:", c.index()
            for indx in self.adopted_index:
                #print ".... indx", indx
                if (indx == rindx):
                    if c.index() == indx:
                        c_tmp       = copy.deepcopy( c )
                        c_tmp.owner = self.rank

                        print "   ...adopting ({},{})".format(c_tmp.index()[0],c_tmp.index()[1])
                        self.cells  = cappend( self.cells, c_tmp )

        for i, c in enumerate( self.virtuals ):
            if c.index() == rindx:
                #self.virtuals = np.delete( self.virtuals, i )
                self.virtuals = cdel( self.virtuals, i )
                break


        #self.adopted_index = []
        #self.adopted_parent = []



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


def cappend(arr, x):
    #arr.append(x)
    #return arr

    tmpa = copy.deepcopy( arr )
    tmpx = copy.deepcopy( x )

    #tmpa = np.append(tmpa, tmpx)
    #return tmpa

    #return tmpa.append( tmpx )
    #return arr.append(x)

    tmpa.append( tmpx )
    return tmpa

    #return np.append(arr, x)

def cdel(arr, i):
    #del arr[i]
    #return arr

    tmp = copy.deepcopy( arr )

    #tmp = np.delete(tmp, i)
    del tmp[i]
    return tmp

    #del arr[i]
    #return arr

    #del tmp[i]
    #return tmp

    #return np.delete(arr, i)


#General communication routines for inter-cell needs
def communicate(nodes):

    #pack for sending
    for node in nodes:
        node.pack_virtuals()

    for node in nodes:
        node.clear_virtuals()
    
    #receive
    for node in nodes:
        for k, c in enumerate(node.send_queue_cells):
            for o in node.send_queue_address[k]:
                #print "sending {} to {}".format(k, o)
                nodes[o].virtuals = cappend( nodes[o].virtuals, c )

    #clear queue lists
    for node in nodes:
        node.clear_queue()


def adopt(nodes):

    #rank virtuals for adoption
    #for node in nodes:
    #    node.rank_virtuals()

    #randomly pick who gets to adopt
    nodes[ np.random.randint(Nrank) ].rank_virtuals()
    #nodes[0].rank_virtuals()


    #communicate adoptions and judge if they are good
    adopted = []
    for node in nodes:
        print node.rank, " has purgelist of ", node.purged
        for (indx, parent) in zip( node.adopted_index, node.adopted_parent):
            if not(indx in adopted):
                print "{} adopting ({},{}) from {}".format(node.rank, indx[1], indx[0], parent)

                node.adopt(indx)

                #print "     0", nodes[0].purged
                #print "     1", nodes[1].purged
                #print "     2", nodes[2].purged

                #nodes[parent].purged.append( indx )
                nodes[parent].purged = cappend(nodes[parent].purged, indx)

                #print "   parent = ", parent, " purgelist:", nodes[parent].purged
                #print "    +0", nodes[0].purged
                #print "    +1", nodes[1].purged
                #print "    +2", nodes[2].purged

                grid[indx] = node.rank
                adopted.append( indx )
        #print node.rank, " after has purgelist of ", node.purged

    #print " "

    #purge adopted children
    for node in nodes:
        print "purging...", node.rank, " cells:", node.purged
        node.purge()






##################################################    
# setup configuration
xmin = ymin = 0.0
xmax = ymax = 1.0

Nrank = 3
Nx = 10
Ny = 10

grid = np.zeros( (Nx, Ny) )




##################################################
# main program loop
def main():

    ## Plot and figure set-up
    fig = figure(figsize=(8,8))
    rc('font', family='serif', size=12)
    rc('xtick')
    rc('ytick')
    
    gs = GridSpec(2, 2)
    gs.update(hspace = 0.5)
    
    
    axs = []
    axs.append( subplot(gs[0]) )
    axs.append( subplot(gs[1]) )
    axs.append( subplot(gs[2]) )
    axs.append( subplot(gs[3]) )

    
    ##################################################
    # Initialize ranks (this is done behind the scenes 
    # like in real MPI startup
    np.random.seed(3)

    nodes = []
    #nodes = np.array([], dtype=np.object)

    for i in range(Nx):
        for j in range(Ny):
            val = np.random.randint(Nrank)
            grid[i,j] = np.float(val)

    # load nodes
    for rank in range(Nrank):
        n = node(rank)
    
        for i in range(Nx):
            for j in range(Ny):
                if grid[i,j] == rank:
                    c = cell(i, j, rank)
                    #n.cells = np.append(n.cells, c)
                    n.cells = cappend(n.cells, c)
                    #n.cells.append( c )
        nodes.append(n)
        #nodes = cappend(nodes, n)
        #nodes = np.append(nodes, n)

    #communicate(nodes)
    ################################################## 


    for t in range(100):
    #for t in [0]:
    #for t in [0,1]:
        print "t=",t

        for s in range(10):
            communicate(nodes)
            adopt(nodes)

        plot_grid(axs[0], nodes)
        plot_node(axs[1], nodes[0])
        plot_node(axs[2], nodes[1])
        plot_node(axs[3], nodes[2])

        pause(1.0)



    show()
    #draw()
    #pause(10.0)





if __name__ == "__main__":
    main()
