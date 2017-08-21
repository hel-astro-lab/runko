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
        #label = str(Nv)
        label = "({},{})/{}".format(i,j,Nv)
        ax.text(jy, ix, label, ha='center',va='center')


    for c in n.virtuals:
        (i, j) = c.index()
        ix = conf.xmin + conf.xmax*(i+0.5)/conf.Nx
        jy = conf.ymin + conf.ymax*(j+0.5)/conf.Ny
        label = "Vir"
        ax.text(jy, ix, label, ha='center',va='center')

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



# open memory slot for put/get
Nqueue = 200
virtual_indices = -1*np.ones( (Nqueue, 4), dtype=int) #(id, i, j, owner)
win_virtual_indices = MPI.Win.Create( virtual_indices, comm=comm)



def transfer_cells(cells, indexes, dest):

    #print rank, " transferring..."

    recv_buffer = -1*np.ones( (Nqueue, 4), dtype=int)
    for k, indx in enumerate(indexes):
        c = cells[indx]
        (i,j) = c.index()
        owner = c.owner
        recv_buffer[k, :] = [rank, i, j, owner]

    win_virtual_indices.Put( [recv_buffer, MPI.INT], dest)


def unpack_incoming_cells():
    k = 0
    for (cid, i, j, owner) in virtual_indices:
        if cid != -1:
            #print "{}:   {}, ({},{}) {}".format(rank, cid, i, j, owner)
            c = cell(i, j, owner)
            n.virtuals = cappend( n.virtuals, c)

            virtual_indices[k, :] = -1 #clean the list
        k += 1
    

def communicate():

    for myid in range(Nrank):
    #for myid in [1]:
        #print "-------sender is {}".format(myid)
        #print myid, " queue:", n.send_queue_address

        #if rank == 1:
        #    for i, c in enumerate(n.send_queue_cells):
        #        #print "in queue: ", c.index(), " for ",n.send_queue_address[i]

        win_virtual_indices.Fence()

        #send
        if myid == rank:
            for dest in range(Nrank):
            #for dest in [0]:
                if dest == rank:
                    continue

                #print "-------to : {}".format(dest)
                indx = []
                for i, address in enumerate(n.send_queue_address):
                    #print i, "address", address
                    if dest in address:
                        indx.append( i ) 
                #print indx
                transfer_cells( n.send_queue_cells, indx, dest )


        #unpack and clean afterwards
        win_virtual_indices.Fence()
        unpack_incoming_cells()
    n.clear_queue()



##################################################
# first communication
plot_node(axs[0], n, 0)
n.pack_all_virtuals()  
n.clear_virtuals()
communicate()



#win_virtual_indices.Fence()
#if master:
#
#    print n.send_queue_address
#
#    for dest in range(Nrank):
#        if dest == rank:
#            continue
#
#        indx = []
#        for i, address in enumerate(n.send_queue_address):
#            if dest in address:
#                indx.append( i ) 
#        print "indices to be sent: ", indx
#        print "len:", len(indx)
#
#    #print "sending: ", n.send_queue_cells[0].index()
#    #print " to :", n.send_queue_address[0]
#        
#        transfer_cells( n.send_queue_cells, indx, dest )
#
#    #for dest in n.send_queue_address[0]:
#    #    print "    to->", dest
#    #    transfer_cell( n.send_queue_cells[0], dest)
#win_virtual_indices.Fence()
#
#print "{}: virtual_indices is: ({},{},{})".format(rank, 
#        virtual_indices[0,1], 
#        virtual_indices[0,2], 
#        virtual_indices[0,3])
#unpack_incoming_cells()




plot_node(axs[0], n, 1)

# initial load balance
#for t in range(10):
#    plot_node(axs[0], n, t)
#
#    n.pack_virtuals()  
#    n.clear_virtuals()



# now start simulation




win_virtual_indices.Free()














