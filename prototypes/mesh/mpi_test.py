import numpy as np
import os, sys
from mpi4py import MPI

import conf
from logi import cappend, cdel 
from logi import cell, node


##################################################    
# MPI start-up
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
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
plt.fig = plt.figure(1, figsize=(8,8))
plt.rc('font', family='serif', size=12)
plt.rc('xtick')
plt.rc('ytick')

gs = plt.GridSpec(1, 1)
gs.update(hspace = 0.5)

axs = []
axs.append( plt.subplot(gs[0]) )




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
              #vmax = Nrank,
              #alpha=0.5
              )

def plot_node(ax, n, lap):
    tmp_grid = np.ones( (conf.Nx, conf.Ny) ) * -1.0

    for c in n.cells:
        (i, j) = c.index()
        tmp_grid[i,j] = c.owner

    for c in n.virtuals:
        (i,j) = c.index()
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
    fname = conf.path + '/node_{}_{}.png'.format(rank, lap)
    plt.savefig(fname)




##################################################    
# mpiGrid setup and configuration
conf.xmin = conf.ymin = 0.0
conf.xmax = conf.ymax = 1.0

conf.Nx = 10
conf.Ny = 10

n = node(rank, conf.Nx, conf.Ny)
np.random.seed(4)

# Compute (initial) total work to balance and my personal effort needed
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
            n.cells = cappend(n.cells, c)



# open memory slot for put/get for cell transfer
Nqueue = 200
Ncellp = 7
# holds information of (id, i, j, owner, Nvir, Ncom)
virtual_indices = np.zeros( (Nqueue, Ncellp), dtype=int) 
win_virtual_indices = MPI.Win.Create( virtual_indices, comm=comm)



def transfer_cells(indexes, dest):

    #print rank, " transferring..."

    #fundamental cell building blocks
    #send_buffer = np.zeros( (Nqueue, Ncellp), dtype=int)
    #send_buffer[:,0] = -1

    virtual_indices[:,0] = -1

    for k, indx in enumerate(indexes):
        c = n.send_queue_cells[indx]
        (i,j)  = c.index()
        owner  = c.owner
        Nvir   = c.number_of_virtual_neighbors
        Ncom   = c.communications
        topown = c.top_owner
        virtual_indices[k, :] = [rank, i, j, owner, Nvir, Ncom, topown]
        #send_buffer[k, :] = [rank, i, j, owner, Nvir, Ncom, topown]
    #win_virtual_indices.Put( [send_buffer, MPI.INT], dest)
    win_virtual_indices.Put( [virtual_indices, MPI.INT], dest)
    #comm.Isend(send_buffer, dest=dest)




def unpack_incoming_cells():
    k = 0
    for (cid, i, j, owner, Nvir, Ncom, topown) in virtual_indices:
        #print "reading: ", virtual_indices[k,:]
        if cid != -1:
            print "reading: ", virtual_indices[k,:]

            #check if incoming virtual is actually our own cell
            # this means that adoption has occured and we need to 
            # change this cell into virtual
            to_be_purged = False
            for q, c in enumerate(n.cells):
                if (i,j) == c.index():
                    print "cell ({},{}) seems to be kidnapped from me {}!".format(i,j, rank)

                    to_be_purged = True


            if to_be_purged:
                n.virtuals = cappend( n.virtuals, n.cells[q] )
                n.cells = cdel( n.cells, q ) #remove from real cells

            else:
                #print "{}:   {}, ({},{}) {}".format(rank, cid, i, j, owner)

                c                             = cell(i, j, owner)
                c.number_of_virtual_neighbors = Nvir
                c.communications              = Ncom
                c.top_owner                   = topown

                n.virtuals = cappend( n.virtuals, c)

            virtual_indices[k, 0] = -1 #clean the list
        k += 1
    

# Routine to communicate cells between nodes.
# For now, a simple loop over all senders to all destinations
# is done with WIN.Fence synchronizing the calls.
def communicate():
    for myid in range(Nrank):
        win_virtual_indices.Fence()

        #send
        if myid == rank:
            for dest in range(Nrank):
                if dest == rank:
                    continue

                indx = []
                for i, address in enumerate(n.send_queue_address):
                    #print i, "address", address
                    if dest in address:
                        indx.append( i ) 
                transfer_cells( indx, dest )

        #unpack and clean afterwards
        win_virtual_indices.Fence()
        unpack_incoming_cells()
    n.clear_queue()



# Decide who to adopt
# Every node has their own council, although they should come
# to same conclusion (hopefully)
def adoption_council():

    # count how many I can adopt; based on my 
    # occupance AND general work level of the society
    print "Current workload: {} / ideal: {}".format( len(n.cells), n.work_goal)
    quota = n.work_goal - len(n.cells)
    print "  => adoption quota: {}".format(quota)

    #if we are overworking, do not allow adoption
    if quota < 0.0:
        return


    #quota is positive, now lets adopt!
    Nvirs, Ncoms, Tvirs, owners, index_list = n.rank_virtuals()
    indx_sorted = np.argsort(Nvirs).tolist()

    #pick virtual cell that has largest number of
    # virtual neighbors AND shares values mostly with me

    #adopted_indices = []
    for i in list(reversed(indx_sorted)):
        if Tvirs[i] == rank:
            n.adopted_index.append( index_list[i] )
            n.adopted_parent.append( index_list[i] )

            break
        #print "virtual cell ({},{}): Nvirs: {} | Ncoms: {} | TopO: {} (parent: {})".format(
        #       index_list[i][1],
        #       index_list[i][0],
        #       Nvirs[i], 
        #       Ncoms[i],
        #       Tvirs[i],
        #       owners[i]
        #       )

    print "winner is:"
    print "virtual cell ({},{}): Nvirs: {} | Ncoms: {} | TopO: {} (parent: {})".format(
               index_list[i][0],
               index_list[i][1],
               Nvirs[i], 
               Ncoms[i],
               Tvirs[i],
               owners[i]
               )


    #return adopted_indices, 
    #inform parents that their child is taken; sorry!




##################################################
# first communication
plot_node(axs[0], n, 0)
n.pack_all_virtuals() 
n.clear_virtuals()
communicate()
plot_node(axs[0], n, 1)


#second round with adoption
#adoption_council()
#n.adopt()
#n.pack_all_virtuals() 
#n.clear_virtuals()
#communicate()
#
#
#plot_node(axs[0], n, 2)

#print n.mpiGrid


# initial load balance
#for t in range(1,2):
#    plot_node(axs[0], n, t)
#
#
#
#    n.pack_virtuals()  
#    n.clear_virtuals()



# now start simulation



# Free MPI processes and memory windows
#win_virtual_indices.Free()














