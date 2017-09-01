import numpy as np
import os, sys
from mpi4py import MPI

import conf
from logi import cappend, cdel 
from logi import cell, node
from logi import Grid


import pickle



##################################################    
# MPI start-up
comm  = MPI.COMM_WORLD
rank  = comm.Get_rank()
Nrank = comm.Get_size()

master_rank = 0
master = (rank == master_rank) #make a master flag

print "Hi from rank {}".format(rank)
if master:
    print "master is {}".format(rank)

    # Path to be created (avoid clutter by issuing this with master)
    if not os.path.exists(conf.path):
        os.makedirs(conf.path)



##################################################
#set up total view of the system
import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm
cmap = pal.wesanderson.Moonrise1_5.mpl_colormap



# set up plotting and figure
plt.fig = plt.figure(1, figsize=(8,4))
plt.rc('font', family='serif', size=12)
plt.rc('xtick')
plt.rc('ytick')

gs = plt.GridSpec(2, 1)
gs.update(hspace = 0.5)

axs = []
axs.append( plt.subplot(gs[0]) )
axs.append( plt.subplot(gs[1]) )




# plotting tools
def imshow(ax, grid):
    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(conf.xmin, conf.xmax)
    ax.set_ylim(conf.xmin, conf.xmax)

    extent = [conf.xmin, conf.xmax, conf.ymin, conf.ymax]

    mgrid = np.ma.masked_where(grid == -1.0, grid)
    ax.imshow(mgrid,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = cmap,
              vmin = 0.0,
              vmax = Nrank-1,
              aspect='auto',
              #vmax = Nrank,
              #alpha=0.5
              )

def plot_grid(ax, n, lap):
    imshow(ax, n.mpiGrid)

    #save
    slap = str(lap).rjust(4, '0')
    fname = conf.path + '/grid_{}.png'.format(slap)
    plt.savefig(fname)




def plot_node(ax, n, lap):
    tmp_grid = np.ones( (conf.Nx, conf.Ny) ) * -1.0

    for c in n.cells:
        (i, j) = c.index()
        #check dublicates
        if tmp_grid[i,j] != -1.0:
            print "{}: ERROR in real cells at ({},{})".format(rank, i,j)
            sys.exit()
        tmp_grid[i,j] = c.owner

    for c in n.virtuals:
        (i,j) = c.index()
        if tmp_grid[i,j] != -1.0:
            print "{}: ERROR in virtual cells at ({},{})".format(rank, i,j)
            sys.exit()
            
        tmp_grid[i,j] = c.owner

    imshow(ax, tmp_grid)
    #imshow(ax, n.mpiGrid)

    # add text label about number of neighbors
    for c in n.cells:
        (i, j) = c.index()
        ix = conf.xmin + conf.xmax*(i+0.5)/conf.Nx
        jy = conf.ymin + conf.ymax*(j+0.5)/conf.Ny

        Nv = n.number_of_virtual_neighbors(c)
        label = str(Nv)
        #label = "({},{})/{}".format(i,j,Nv)
        ax.text(jy, ix, label, ha='center',va='center')


    #for c in n.virtuals:
    #    (i, j) = c.index()
    #    ix = conf.xmin + conf.xmax*(i+0.5)/conf.Nx
    #    jy = conf.ymin + conf.ymax*(j+0.5)/conf.Ny
    #    label = "Vir"
    #    ax.text(jy, ix, label, ha='center',va='center')

    ax.set_title(str(len(n.virtuals))+"/"+str(len(n.cells)))


    #save
    slap = str(lap).rjust(4, '0')
    fname = conf.path + '/node_{}_{}.png'.format(rank, slap)
    plt.savefig(fname)



def plot_node_data(ax, n, lap):

    NxFull   = conf.Nx * conf.NxCell
    NyFull   = conf.Ny * conf.NyCell
    fullMesh = np.zeros( (NxFull, NyFull) )

    
    #loop over everything and collect the data
    for c in n.cells:
        (i, j) = c.index()
        i1 = i * conf.NxCell
        j1 = j * conf.NyCell

        mesh = c.grid.mesh
        for ix in range(conf.NxCell):
            for iy in range(conf.NyCell):
                fullMesh[i1 + ix, j1 + iy] = mesh[ix, iy]



    #imshow(ax, fullMesh)
    ax.clear()
    ax.minorticks_on()
    ax.set_xlim(conf.xmin, conf.xmax)
    ax.set_ylim(conf.xmin, conf.xmax)

    extent = [conf.xmin, conf.xmax, conf.ymin, conf.ymax]

    mgrid = np.ma.masked_where(fullMesh == 0.0, fullMesh)
    ax.imshow(mgrid,
              extent=extent,
              origin='lower',
              interpolation='nearest',
              cmap = cm.get_cmap('plasma'),
              vmin = 0.0,
              vmax = 1.0,
              aspect='auto',
              #alpha=0.5
              )
    #save
    slap = str(lap).rjust(4, '0')
    #fname = conf.path + '/data_{}_{}.png'.format(rank, slap)
    #plt.savefig(fname)





##################################################    
# mpiGrid setup and configuration
conf.xmin = conf.ymin = 0.0
conf.xmax = 1.0
conf.ymax = 1.0

conf.Nx = 10
conf.Ny = 40


conf.NxCell = 2
conf.NyCell = 2


n = node(rank, conf.Nx, conf.Ny)
np.random.seed(4)

# Compute (initial) total work in order to balance and my personal effort needed
n.total_work = conf.Nx * conf.Ny
n.work_goal  = n.total_work / Nrank


if master:
    print "Average work load is {} cells of {}".format(n.work_goal, n.total_work)



# initialize grid ownership (master does this, for now)
if master:
    for i in range(conf.Nx):
        for j in range(conf.Ny):
            val = np.random.randint(Nrank)
            n.mpiGrid[i,j] = np.float(val)
n.mpiGrid = comm.bcast( n.mpiGrid, root=master_rank)
    

#load cells into local node
for i in range(conf.Nx):
    for j in range(conf.Ny):
        if n.mpiGrid[i,j] == rank:
            c = cell(i, j, rank)

            #load data into cell
            grid = Grid(0.0, 1.0, 0.0, 1.0, conf.NxCell, conf.NyCell)
            c.grid = grid


            n.cells = cappend(n.cells, c)


#from subgrid index (i,j) to global coordinates
def ij2xy(i,j):
    dx = (conf.xmax - conf.xmin)/(conf.Nx * conf.NxCell)
    dy = (conf.ymax - conf.ymin)/(conf.Ny * conf.NyCell)

    x = conf.xmin + i*dx
    y = conf.ymin + j*dy

    return x,y

def gauss(x,y):
    g = np.exp(-((y - 0.25)/0.2 )**2)

    return g

    

#insert some data
for c in n.cells:
    (i,j) = c.index()

    for ii in range(conf.NxCell):
        for jj in range(conf.NyCell):
            yi, xi = ij2xy(i+ii, j+jj) #we flip the coordinates here

            c.grid.mesh[ii, jj] = xi
            #c.grid.mesh[ii, jj] = gauss(xi, yi)
    




def send_virtual_cells_pickle(indexes, dest):
    #use normal python array and Cell class 
    # pickle will turn this into bytestream for sending
    for k, indx in enumerate(indexes):
        comm.isend( n.send_queue_cells[indx], dest=dest, tag=data_tag)


# Purge here actually means shifting from real to virtual state
def purge_kidnapped_cells():

    for (indx,kidnapper) in zip( n.kidnap_index, n.kidnapper ):
        #print "inside purge... ",rank,indx
        for q, c in enumerate(n.cells):
            if c.index() == indx:
                print "{}: cell ({},{}) seems to be kidnapped from me by {}!".format(rank, indx[0], indx[1], kidnapper)

                n.virtuals = cappend( n.virtuals, c)
                n.mpiGrid[c.i, c.j] = kidnapper

                #and remove from; 
                n.cells = cdel( n.cells, q ) #remove from real cells

                break



def unpack_incoming_cell(inc_cell):

    for q, c in enumerate(n.cells):
        if inc_cell.index() == c.index():
            print "{}: ERROR: I already have cell ({},{}) from {}!".format(rank, inc_cell.i, inc_cell.j, inc_cell.owner)
            #XXX this is where we could compare cells and pick the newest
            
    #this should be update, not replace
    n.virtuals = cappend( n.virtuals, inc_cell)

   
comm_tag  = 0
data_tag  = 1
adopt_tag = 2


# Routine to communicate cells between nodes.
# For now, a simple loop over all senders to all destinations
def communicate_send_cells():

    for dest in range(Nrank):
        if dest == rank:
            continue

        indx = []
        for i, address in enumerate(n.send_queue_address):
            #print i, "address", address
            if dest in address:
                indx.append( i ) 

        #print "{}: sending {} cells to {} ".format(rank, len(indx), dest)


        #initial message informing how many cells are coming
        Nincoming_cells = len(indx)
        comm.isend(Nincoming_cells, dest=dest, tag=comm_tag)


        #send_virtual_cells( indx, dest ) #numpy
        send_virtual_cells_pickle( indx, dest ) #pickle


    #n.clear_queue()
    return


def communicate_recv_cells():
    for source in range(Nrank):
        if source == rank:
            continue

        #Communicate with how many cells there are incoming
        req = comm.irecv(source=source, tag=comm_tag)
        Nincoming_cells = req.wait()
        #print " First contact: we are expecting {} cells from {}".format(Nincoming_cells, source)

        # create a list of pending receives
        reqs = [] 
        for irecvs in range(Nincoming_cells):
            req = comm.irecv(source=source, tag=data_tag)
            reqs.append( req )


        #process the queue
        Nrecs = 0
        while Nrecs < Nincoming_cells:

            for irecvs in range( len(reqs) ):
                req = reqs[ irecvs ]
                inc_cell = req.wait() #XXX remove this barrier

                #print "   received 1 cell ({},{})".format(inc_cell.i, inc_cell.j)
                unpack_incoming_cell( inc_cell )

                Nrecs += 1

        #print "{} successfully received and unpacked everything".format(rank)


    return 


#blocking version of communicate
def communicate_cells():
    communicate_send_cells()
    communicate_recv_cells()


# Decide who to adopt
# Every node has their own council, although they should come
# to same conclusion (hopefully)
def adoption_council():

    # count how many I can adopt; based on my 
    # occupance AND general work level of the society
    print "{}: Current workload: {} / ideal: {}".format( rank, len(n.cells), n.work_goal)
    quota = n.work_goal - len(n.cells)
    print "  => {}: adoption quota: {}".format(rank, quota)

    #if we are overworking, do not allow adoption
    #if quota < 0:
    #    return

    quota = np.clip(quota, 2, None)


    #quota is positive, now lets adopt!
    Nvirs, Ncoms, Tvirs, owners, index_list = n.rank_virtuals()
    indx_sorted = np.argsort(Nvirs).tolist()

    #pick virtual cell that has largest number of
    # virtual neighbors AND shares values mostly with me

    #limit the quota
    #if quota > 3:
    #    quota = 3
    #quota = 1


    Nadopts = 0
    for i in list(reversed(indx_sorted)):
        if Tvirs[i] == rank:
            n.adopted_index.append( index_list[i] )
            n.adopted_parent.append( owners[i] )
            
            Nadopts += 1

            print "{}: got adoption command for ({},{}) of {}".format(rank, index_list[i][0], index_list[i][1], owners[i])


        print "virtual cell ({},{}): Nvirs: {} | Ncoms: {} | TopO: {} (parent: {})".format(
               index_list[i][0],
               index_list[i][1],
               Nvirs[i], 
               Ncoms[i],
               Tvirs[i],
               owners[i]
               )

        if Nadopts >= quota:
            break

    return


def communicate_send_adoptions():
    for dest in range(Nrank):
        if dest == rank:
            continue

        kidnaps = []
        for q, address in enumerate( n.adopted_parent ):
            #if address == dest:
            print "{}: sending info about my kidnap to {}".format(rank, dest)
            kidnaps.extend( (n.adopted_index[q]) )



        print "sending this: ", kidnaps
        comm.isend(kidnaps, dest=dest, tag=adopt_tag)


def communicate_recv_adoptions():
    n.kidnap_index = []
    n.kidnapper    = []
    for source in range(Nrank):
        if source == rank:
            continue
        
        adopt_req = comm.irecv(source=source, tag=adopt_tag)
        kidnaps = adopt_req.wait()
        print "receiving this:", kidnaps

        for i in range(0,len(kidnaps), 2):
            indx = ( kidnaps[i], kidnaps[i+1] )  #tuplify

            #check if it is mine
            if n.is_local( indx ):
                print "{}: oh gosh, my child ({},{}), has been kidnapped by evil {}".format(rank, indx[0], indx[1], source)
                n.kidnap_index.append( indx )
                n.kidnapper.append( source )
            else:
                print "{}: BUU!, {} kidnapped ({},{})".format(rank, source, indx[0], indx[1])
                n.mpiGrid[ indx[0], indx[1] ] = source



#blocking version of adoption communicate
def communicate_adoptions():
    communicate_send_adoptions()
    communicate_recv_adoptions()





##################################################
# first communication
plot_node_data(axs[1], n, 0)
plot_node(axs[0], n, 0)

n.pack_all_virtuals() 
n.clear_virtuals()
communicate_send_cells()
communicate_recv_cells()


#plot_node(axs[0], n, 1)
if master:
    plot_grid(axs[0], n, 1)


#second round with adoption
adoption_council()
n.adopt()
communicate_adoptions()
purge_kidnapped_cells()
n.pack_all_virtuals() 
n.clear_virtuals()
communicate_cells()



#plot_node(axs[0], n, 2)
if master:
    plot_grid(axs[0], n, 2)


#sys.exit()

# initial load balance burn-in
for t in range(3,25):

    adoption_council()
    n.adopt()
    communicate_adoptions()
    purge_kidnapped_cells()

    n.pack_all_virtuals() 
    n.clear_virtuals()

    communicate_cells()

    plot_node(axs[0], n, t)
    plot_node_data(axs[1], n, 1)

    if master:
        plot_grid(axs[0], n, t)




# now start simulation

















