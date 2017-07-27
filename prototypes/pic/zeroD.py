import numpy as np
import math
from pylab import *
import os, sys
from scipy.stats import gaussian_kde
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import kn

#set seed to get reproducible errors & results
np.random.seed(1)

#update system path as we are on the subdirectory
sys.path.append('/Users/natj/projects/radpic/prototype')
import radpic as rpic


#dimensionless inputs
##################################################
rpic.Nppc = 200 #particles per cell
rpic.delgamma = 0.5 #(k T/m c^2, T electron/ion temperature

rpic.qe = 1.0 #electron charge
rpic.qi = -1.0 #positron 

rpic.me = 1.0 #electron mass
rpic.mi = 1.0 #ion mass (or actually mass to charge ratio)
rpic.gamma = 0.0 #flow drift gamma (gamma0) #now sets beta = v/c  because < 1
rpic.delta = 10.0 #plasma skin depth
rpic.Te_Ti = 1.0 #T_e / T_i
rpic.c = 1.0 #Computational speed of light


#two or three D
rpic.twoD = True


#grid dimensions
rpic.Nx = 2
rpic.Ny = 2
rpic.Nz = 1


rpic.Nx_wrap = True
rpic.Ny_wrap = True
rpic.Nz_wrap = True

rpic.Np = rpic.Nx * rpic.Nppc #total number of particles


#initialize according to skin depth/cell
rpic.grid_scale_skindepth(rpic.delta)


#Initialize the grid
rpic.init()


#Mesh grid for plotting
#XX, YY, ZZ = np.meshgrid(np.linspace(rpic.grid_xmin, rpic.grid_xmax, rpic.Nx),
#                         np.linspace(rpic.grid_ymin, rpic.grid_ymax, rpic.Ny),
#                         np.linspace(rpic.grid_zmin, rpic.grid_zmax, rpic.Nz)
#                        )
XX, YY, ZZ = np.meshgrid(np.linspace(0.0, rpic.Nx-1, rpic.Nx),
                         np.linspace(0.0, rpic.Nx-1, rpic.Ny),
                         np.linspace(0.0, rpic.Nx-1, rpic.Nz)
                        )



##################################################
# Path to be created
path = "out"
if not os.path.exists(path):
    os.makedirs(path)



##################################################
#set up figure
fig = figure(figsize=(10, 12), dpi=200)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(3, 3)
gs.update(hspace = 0.2)
gs.update(wspace = 0.2)

ax1 = subplot(gs[0,0])
ax2 = subplot(gs[0,1])
ax3 = subplot(gs[0,2])

ax4 = subplot(gs[1,0])
ax5 = subplot(gs[1,1])
ax6 = subplot(gs[1,2])

ax7 = fig.add_subplot(gs[2,0], projection='3d')
ax8 = fig.add_subplot(gs[2,1], projection='3d')
ax9 = fig.add_subplot(gs[2,2], projection='3d')

#for ax in [ax1, ax2, ax3]:
#    ax.set_xlim(-5, 5)


##################################################
#inject particles

#initialize particles
delgamma_i = rpic.delgamma
delgamma_e = rpic.delgamma *(rpic.mi/rpic.me)*rpic.Te_Ti

for i in range(rpic.Nx):
    for j in range(rpic.Ny):
        for k in range(rpic.Nz):
            cell = rpic.CellClass(i,j,k) 

            #cell limits
            xmin,xmax, ymin,ymax, zmin,zmax = rpic.grid_limits(i,j,k)
            
            xarr = np.linspace(xmin, xmax, rpic.Nppc/2)
            for n in range(rpic.Nppc/2):

                #random uniform location inside cell
                x = xmin + (xmax-xmin)*np.random.rand()
                y = ymin + (ymax-ymin)*np.random.rand()
                z = zmin + (zmax-zmin)*np.random.rand()

                vx, vy, vz, u = rpic.boosted_maxwellian(delgamma_e, rpic.gamma, dims=3)

                w = rpic.qe * 1.0

                cell.particles = np.concatenate((cell.particles, [[x, y, z, vx, vy, vz, w]]), axis=0) #add 6D phase space 


                #Add ion/positron to exact same place
                vxi, vyi, vzi, ui = rpic.boosted_maxwellian(delgamma_i, rpic.gamma, dims=3)

                wi = rpic.qi*1.0


                cell.particles = np.concatenate((cell.particles, [[x, y, z, vxi, vyi, vzi, wi]]), axis=0) #add 6D phase space 


                cell.Npe += 2
            rpic.mpiGrid[i,j,k] = cell





def kde_scipy(x, x_grid, bandwidth='scott', **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    #kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)

    kde = gaussian_kde(x, bw_method=bandwidth, **kwargs)
    return kde.evaluate(x_grid)


def plot_pdf(ax, sample, xlabel='', color='blue'):
    edges = np.linspace(-7, 7, 100)
    
    ax.cla()
    ax.set_xlabel(xlabel)
    ax.set_xlim(-5, 5)

    #kernel density estimator
    pdf1 = kde_scipy(sample, edges, bandwidth='scott')
    ax.plot(edges, pdf1, color=color, alpha=0.8)

    return ax


def plot_pdf2d(ax, xs, ys, xlabel='', ylabel=''):

    ax.cla()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

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


def plot_pdf3d(ax, xs, ys, xlabel='', ylabel=''):

    ax.cla()
    ax.set_ylabel(xlabel)
    ax.set_xlabel(ylabel)
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


def plot_pdf3d_isosurf(ax, xs, ys):

    ax.cla()
    ax.set_ylabel(r'$v_x$')
    ax.set_xlabel(r'$v_y$')
    ax.set_zlabel(r'$v_z$')

    xmin = -6.0
    xmax = 6.0
    ymin = -6.0
    ymax = 6.0

    zmin = -6.0
    zmax = 6.0

    kde = gaussian_kde(values, bw_method='scott')
    xi, yi, zi = np.mgrid[xmin:xmax:50j, ymin:ymax:50j, zmin:zmax:50j]

    # Evaluate the KDE on a regular grid...
    coords = np.vstack([item.ravel() for item in [xi, yi, zi]])
    density = kde(coords).reshape(xi.shape)

    # Visualize the density estimate as isosurfaces
    mlab.contour3d(xi, yi, zi, density, opacity=0.5)
    mlab.axes()
    mlab.show()


    return ax



##################################################
##################################################
##################################################

#initialize fields and currents

#x-dir electric field to accelerate particles
rpic.Ex[0,0,0] = 0.0
rpic.Ey[0,0,0] = 0.0
rpic.Ez[0,0,0] = 0.0

#Non-zero y-dir B-field
rpic.Bx[0,0,0] = 0.0
rpic.By[0,0,0] = 0.0
rpic.Bz[0,0,0] = 0.0


max_steps = 10
for step in range(1, max_steps):
    print " step: ",step


    rpic.Vay_update_velocities(rpic.mpiGrid)

    #do not update E or J
    #rpic.push_E()
    #rpic.conserving_deposit_current(rpic.mpiGrid)

    rpic.sort_particles_between_cells(rpic.mpiGrid)


    #apply filters
    for sweeps in range(32):
        rpic.filter_current(0.5,1) #x sweep
        rpic.filter_current(0.5,2) #y sweep
    #rpic.filter_current(-1.0/6.0, 1) #put some power back with negative sweeps
    #rpic.filter_current(-1.0/6.0, 2)


    #I/O
    ################################################## 

    #collect particles and dissect into elements
    particles = rpic.collect_grid(rpic.mpiGrid)
    electrons, positrons = rpic.divide_species(particles)

    x = electrons[:,0]
    y = electrons[:,1]
    z = electrons[:,2]
    vx = electrons[:,3]
    vy = electrons[:,4]
    vz = electrons[:,5]

    ax1 = plot_pdf(ax1, vx, r'$v_x$')
    ax2 = plot_pdf(ax2, vy, r'$v_y$')
    ax3 = plot_pdf(ax3, vz, r'$v_z$')

    ax4 = plot_pdf2d(ax4, vx, vy, r'$v_x$', r'$v_y$')
    #ax5 = plot_pdf2d(ax5, vy, vz, r'$v_y$', r'$v_z$')
    ax6 = plot_pdf2d(ax6, vx, vz, r'$v_x$', r'$v_z$')

    ax7 = plot_pdf3d(ax7, vx, vy, r'$v_x$', r'$v_y$')
    #ax8 = plot_pdf3d(ax8, vy, vz, r'$v_y$', r'$v_z$')
    ax9 = plot_pdf3d(ax9, vx, vz, r'$v_x$', r'$v_z$')

    #ax6 = plot_pdf2d(ax6, vx, vy)
    #ax6 = plot_pdf3d(ax6, vx, vy)
    #ax6 = plot_pdf3d_isosurf(ax6, vx, vy)

    fname = path+'/zeroD_'+str(step)+'.png'
    savefig(fname)

    #clean figure
    #res1.pop(0).remove()
    #res2.pop(0).remove()
    #res3.pop(0).remove()
    #res4.pop(0).remove()

#end of loop



