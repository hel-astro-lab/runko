import numpy as np
import math
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import os, sys

import radpic as rpic


#physical parameters
rpic.e = 1.0
rpic.me = 1.0
rpic.c = 1.0

rpic.threeD = True


#grid dimensions
rpic.Nx = 5
rpic.Ny = 5
rpic.Nz = 5

rpic.Nx_wrap = True
rpic.Ny_wrap = True
rpic.Nz_wrap = True

rpic.grid_xmin=0.0
rpic.grid_xmax=10.0

rpic.grid_ymin=0.0
rpic.grid_ymax=10.0

rpic.grid_zmin=0.0
rpic.grid_zmax=10.0

rpic.Np = 1000 #total number of particles


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
#fig = figure(figsize=(10, 8), dpi=200)
fig = plt.figure()
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

ax1 = fig.add_subplot(224, projection='3d')
ax2 = fig.add_subplot(221, projection='3d')
ax3 = fig.add_subplot(222, projection='3d')
ax4 = fig.add_subplot(223, projection='3d')

ax2.set_title(r'$\bar{E}$')
ax3.set_title(r'$\bar{B}$')
ax4.set_title(r'Current $\bar{J}$')

for ax in [ax1,ax2,ax3,ax4]:
    ax.set_xlim((rpic.grid_xmin, rpic.grid_xmax))
    ax.set_ylim((rpic.grid_ymin, rpic.grid_ymax))
    ax.set_zlim((rpic.grid_zmin, rpic.grid_zmax))
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$z$')




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

                cell.electrons = np.concatenate((cell.electrons, [[x, y, z, vx, vy, vz]]), axis=0) #add 6D phase space 

                cell.Npe += 1

            rpic.mpiGrid[i,j,k] = cell



##################################################
##################################################
##################################################

#initialize fields and currents
rpic.deposit_current(rpic.mpiGrid)
rpic.Yee_currents()


max_steps = 10
for step in range(1, max_steps):
    print " step: ",step

    rpic.push_half_B()

    rpic.update_velocities(rpic.mpiGrid)
    
    rpic.sort_particles_between_cells(rpic.mpiGrid)

    rpic.push_half_B()

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

    res1  = ax1.plot(x, y, z, "k.", markersize=0.5)

    res2 =  ax2.quiver(XX, YY, ZZ, rpic.Ex, rpic.Ey, rpic.Ez, pivot='tail')
    res3 =  ax3.quiver(XX, YY, ZZ, rpic.Bx, rpic.By, rpic.Bz, pivot='tail')
    res4 =  ax4.quiver(XX, YY, ZZ, rpic.Jx, rpic.Jy, rpic.Jz, pivot='tail')



    fname = path+'/threeD_'+str(step)+'.png'
    savefig(fname)

    #clear figure
    res1.pop(0).remove()
    ax2.collections.remove(res2)
    ax3.collections.remove(res3)
    ax4.collections.remove(res4)

#end of loop



