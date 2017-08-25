import numpy as np
import os, sys
from mpi4py import MPI

import conf
from logi import cappend, cdel 
from logi import cell, node


import pickle



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
    slap = str(lap).rjust(4, '0')
    fname = conf.path + '/node_{}_{}.png'.format(rank, slap)
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


def send_virtual_cells_pickle(indexes, dest):
    #use normal python array and Cell class 
    # pickle will turn this into bytestream for sending
    for k, indx in enumerate(indexes):
        comm.isend( n.send_queue_cells[indx], dest=dest, tag=data_tag)



def unpack_incoming_cell(inc_cell):

    #check if incoming virtual is actually our own cell
    # this means that adoption has occured and we need to 
    # change this cell into virtual
    to_be_purged = False
    for q, c in enumerate(n.cells):
        if inc_cell.index() == c.index():
            print "cell ({},{}) seems to be kidnapped from me {}!".format(inc_cell.i, inc_cell.j, rank)

            to_be_purged = True
            break


    if to_be_purged:
        #TODO: choose incoming or current cell?
        n.virtuals = cappend( n.virtuals, inc_cell)
        #n.virtuals = cappend( n.virtuals, n.cells[q] )
        n.mpiGrid[inc_cell.i, inc_cell.j] = inc_cell.owner

        #and remove from; 
        n.cells = cdel( n.cells, q ) #remove from real cells
    else:
        #this should be update, not replace
        n.virtuals = cappend( n.virtuals, inc_cell)

        
   
comm_tag = 0
data_tag = 1


# Routine to communicate cells between nodes.
# For now, a simple loop over all senders to all destinations
def communicate_send():

    for dest in range(Nrank):
        if dest == rank:
            continue

        indx = []
        for i, address in enumerate(n.send_queue_address):
            #print i, "address", address
            if dest in address:
                indx.append( i ) 

        print "{}: sending {} cells to {} ".format(rank, len(indx), dest)


        #initial message informing how many cells are coming
        Nincoming_cells = len(indx)
        comm.isend(Nincoming_cells, dest=dest, tag=comm_tag)


        #send_virtual_cells( indx, dest ) #numpy
        send_virtual_cells_pickle( indx, dest ) #pickle


    #n.clear_queue()
    return


def communicate_receive():
    for source in range(Nrank):
        if source == rank:
            continue

        req = comm.irecv(source=source, tag=comm_tag)
        Nincoming_cells = req.wait()
        print " First contact: we are expecting {} cells from {}".format(Nincoming_cells, source)


        # create a list pending receives
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

        print "{} successfully received and unpacked everything".format(rank)


    return 


#blocking version of communicate
def communicate():
    communicate_send()
    communicate_receive()


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

    #limit the quota
    #if quota > 3:
    #    quota = 3


    Nadopts = 0
    for i in list(reversed(indx_sorted)):
        if Tvirs[i] == rank:
            n.adopted_index.append( index_list[i] )
            n.adopted_parent.append( index_list[i] )
            
            Nadopts += 1

        if Nadopts >= quota:
            break

        #print "virtual cell ({},{}): Nvirs: {} | Ncoms: {} | TopO: {} (parent: {})".format(
        #       index_list[i][1],
        #       index_list[i][0],
        #       Nvirs[i], 
        #       Ncoms[i],
        #       Tvirs[i],
        #       owners[i]
        #       )

    return




##################################################
# first communication
plot_node(axs[0], n, 0)
n.pack_all_virtuals() 
n.clear_virtuals()
communicate_send()
communicate_receive()
plot_node(axs[0], n, 1)


#second round with adoption
adoption_council()
n.adopt()
n.pack_all_virtuals() 
n.clear_virtuals()
communicate()
plot_node(axs[0], n, 2)



# initial load balance burn-in
#for t in range(3,10):
#    adoption_council()
#    n.adopt()
#    n.pack_all_virtuals() 
#    n.clear_virtuals()
#    communicate()
#
#    plot_node(axs[0], n, t)
#    #if master:
#    #    plot_node(axs[0], n, t)




# now start simulation

















