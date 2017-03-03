import numpy as np
import math
from pylab import *
import os, sys

import radpic as rpic


#physical parameters
rpic.e = 1.0
rpic.me = 1.0
rpic.c = 1.0

rpic.twoD = True


#grid dimensions
rpic.Nx = 50
rpic.Ny = 50
rpic.Nz = 1

rpic.Nx_wrap = True
rpic.Ny_wrap = True
rpic.Nz_wrap = True

rpic.grid_xmin=0.0
rpic.grid_xmax=10.0

rpic.grid_ymin=0.0
rpic.grid_ymax=10.0

rpic.grid_zmin=0.0
rpic.grid_zmax=10.0

rpic.Np = 10000 #total number of particles


#Initialize the grid
rpic.init()


#Mesh grid for plotting
XX, YY, ZZ = np.meshgrid(np.linspace(rpic.grid_xmin, rpic.grid_xmax, rpic.Nx),
                         np.linspace(rpic.grid_ymin, rpic.grid_ymax, rpic.Ny),
                         np.linspace(rpic.grid_zmin, rpic.grid_zmax, rpic.Nz)
                        )


##################################################
# Path to be created
path = "out"
if not os.path.exists(path):
    os.makedirs(path)

#set up figure
fig = figure(figsize=(10, 8), dpi=200)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(2, 2)
gs.update(hspace = 0.2)
gs.update(wspace = 0.2)

ax1 = subplot(gs[0,0])

ax2 = subplot(gs[1,0])
ax2.set_title(r'$\bar{E}$')

ax3 = subplot(gs[1,1])
ax3.set_title(r'$\bar{B}$')

ax4 = subplot(gs[0,1])
ax4.set_title(r'Current $J_x$')


for ax in [ax1, ax2, ax3, ax4]:
    ax.set_xlim((rpic.grid_xmin, rpic.grid_xmax))
    ax.set_ylim((rpic.grid_ymin, rpic.grid_ymax))
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')



##################################################
#inject particles
for i in range(rpic.Nx):
    for j in range(rpic.Ny):
        for k in range(rpic.Nz):
            cell = rpic.CellClass() 

            #cell limits
            xmin,xmax, ymin,ymax, zmin,zmax = rpic.grid_limits(i,j,k)
            

            for n in range(rpic.Np/rpic.Ncells):
                #random uniform location inside cell
                x = xmin + (xmax-xmin)*np.random.rand()
                y = ymin + (ymax-ymin)*np.random.rand()
                z = zmin + (zmax-zmin)*np.random.rand()

                #Temperature
                if (n % 2 == 0):
                    vb = 2.0
                else:
                    vb = -2.0

                vx = rpic.Maxwellian(vb)
                vy = rpic.Maxwellian(vb)
                vz = rpic.Maxwellian(vb)

                cell.particles = np.concatenate((cell.particles, [[x, y, z, vx, vy, vz]]), axis=0) #add 6D phase space 

                cell.Npe += 1

            rpic.mpiGrid[i,j,k] = cell



##################################################
##################################################
##################################################

#initialize fields and currents
rpic.deposit_current(rpic.mpiGrid)
rpic.Yee_currents()


max_steps = 3
for step in range(1, max_steps):
    print " step: ",step

    #push_half_B()

    rpic.update_velocities(rpic.mpiGrid)
    
    rpic.sort_particles_between_cells(rpic.mpiGrid)

    #push_half_B()

    rpic.push_E()

    rpic.deposit_current(rpic.mpiGrid)

    rpic.Yee_currents()

    #apply filters

    #I/O
    ################################################## 
    electrons = rpic.collect_grid(rpic.mpiGrid)
    x = electrons[:,0]
    y = electrons[:,1]
    z = electrons[:,2]
    vx = electrons[:,3]
    vy = electrons[:,4]
    vz = electrons[:,5]

    res1 = ax1.plot(x, y, "k.", markersize=0.5)

    #res2 = ax2.plot(XX[0,:,0], rpic.Ex[:,0,0], "b-")
    #res3 = ax3.plot(XX[0,:,0], rpic.By[:,0,0], "b-")
    #res3 = ax3.plot(XX[0,:,0], rpic.Bz[:,0,0], "r-")
    #res4 = ax4.plot(XX[0,:,0], rpic.Jx[:,0,0], "b-")

    fname = path+'/twoD_'+str(step)+'.png'
    savefig(fname)

    #clean figure
    res1.pop(0).remove()
    #res2.pop(0).remove()
    #res3.pop(0).remove()
    #res4.pop(0).remove()

#end of loop



