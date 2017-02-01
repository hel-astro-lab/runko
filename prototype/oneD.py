import numpy as np
import math
from pylab import *
import os, sys
from scipy.stats import gaussian_kde
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import radpic as rpic


#physical parameters
rpic.e = 1.0
rpic.me = 1.0
rpic.c = 1.0


rpic.oneD = True #dimensionality
rpic.twoD = False
rpic.threeD = False


#grid dimensions
rpic.Nx = 100
rpic.Ny = 1
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
fig = figure(figsize=(10, 12), dpi=200)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(3, 2)
gs.update(hspace = 0.2)
gs.update(wspace = 0.2)

ax1 = subplot(gs[0,0])
ax1.set_ylim((-5, 5))
ax1.set_ylabel(r'velocity $v$')

ax2 = subplot(gs[1,0])
ax2.set_ylim((-15, 15))
ax2.set_ylabel(r'$E_x$')

ax3 = subplot(gs[1,1])
ax3.set_ylabel(r'$B_y$, $B_z$')

ax4 = subplot(gs[0,1])
ax4.set_ylim((-5, 5))
ax4.set_ylabel(r'Current $J_x$')

ax5 = subplot(gs[2,0])

ax6 = fig.add_subplot(gs[2,1], projection='3d')
#ax6 = subplot(gs[2,1])



for ax in [ax1, ax2, ax3, ax4]:
    ax.set_xlim((rpic.grid_xmin, rpic.grid_xmax))
    ax.set_xlabel(r'$x$')



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


def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    #kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    #
    # normal scott kernel bandwidth estimator
    kde = gaussian_kde(x, bw_method='scott', **kwargs)
    return kde.evaluate(x_grid)


def slice_pdf(ax, particles):

    N_slices = 10
    cm_subsection = np.linspace(0.0, 1.0, N_slices) 
    colors = [cm.Dark2(x) for x in cm_subsection ]

    ax.cla()
    ax.set_ylabel(r'pdf slices')
    ax.set_xlabel(r'$v_x$')

    
    xlims = np.linspace(rpic.grid_xmin, rpic.grid_xmax, N_slices+1)
    for i in range(len(xlims)-1):
        #print "slice lims", xlims[i], xlims[i+1]
        indx1 = np.where(particles[:,0] > xlims[i])
        indx2 = np.where(particles[:,0] < xlims[i+1])
        indx = np.intersect1d(indx1, indx2)

        vx_sample = particles[indx, 3]
        ax = plot_pdf(ax, vx_sample, color=colors[i])

    return ax


def plot_pdf(ax, sample, color='blue'):
    edges = np.linspace(-7, 7, 100)
    
    #kernel density estimator
    pdf1 = kde_scipy(sample, edges, bandwidth=0.2)
    ax.plot(edges, pdf1, color=color, alpha=0.8)

    return ax


def plot_pdf2d(ax, xs, ys):

    ax.cla()
    ax.set_ylabel(r'$v_x$')
    ax.set_xlabel(r'$v_y$')
    #ax.set_zlabel(r'pdf')

    xmin = -6.0
    xmax = 6.0
    ymin = -6.0
    ymax = 6.0
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([xs, ys])

    kernel = gaussian_kde(values, bw_method='scott')
    #kernel = stats.gaussian_kde(values)

    Z = np.reshape(kernel(positions).T, X.shape)

    #ax.plot_surface(X, Y, Z)
    cfset = ax.contourf(X,Y,Z, cmap='Blues')
    cset = ax.contour(X,Y,Z, colors='k')
    #ax.clabel(cset, inline=1)

    return ax


def plot_pdf3d(ax, xs, ys):

    ax.cla()
    ax.set_ylabel(r'$v_x$')
    ax.set_xlabel(r'$v_y$')
    ax.set_zlabel(r'pdf')

    xmin = -6.0
    xmax = 6.0
    ymin = -6.0
    ymax = 6.0
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([xs, ys])

    kernel = gaussian_kde(values, bw_method='scott')

    Z = np.reshape(kernel(positions).T, X.shape)

    ax.plot_surface(X, Y, Z, cmap='Blues')

    return ax

##################################################
##################################################
##################################################

#initialize fields and currents
rpic.deposit_current(rpic.mpiGrid)
rpic.Yee_currents()


max_steps = 10
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

    res1 = ax1.plot(x, vx, "k.")

    res2 = ax2.plot(XX[0,:,0], rpic.Ex[:,0,0], "b-")
    res3 = ax3.plot(XX[0,:,0], rpic.By[:,0,0], "b-")
    res3 = ax3.plot(XX[0,:,0], rpic.Bz[:,0,0], "r-")
    res4 = ax4.plot(XX[0,:,0], rpic.Jx[:,0,0], "b-")

    ax5 = slice_pdf(ax5, electrons)

    #ax6 = plot_pdf2d(ax6, vx, vy)
    ax6 = plot_pdf3d(ax6, vx, vy)

    fname = path+'/oneD_'+str(step)+'.png'
    savefig(fname)

    #clean figure
    res1.pop(0).remove()
    res2.pop(0).remove()
    res3.pop(0).remove()
    res4.pop(0).remove()


#end of loop



