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
sys.path.append('/Users/natj/projects/plasma/prototypes/pic')
import radpic as rpic


#new dimensionless inputs
##################################################
rpic.Nppc = 64 #particles per cell
#rpic.Nppc = 10 #particles per cell
rpic.delgamma = 2.0e-5 #(k T/m c^2, T electron/ion temperature

rpic.qe = 1.0 #electron charge
rpic.qe = 1.0 #electron

rpic.me = 1.0 #electron mass
rpic.mi = 1.0 #ion mass (or actually mass to charge ratio)
rpic.gamma = 0.5 #flow drift gamma (gamma0) #now sets beta = v/c  because < 1
rpic.delta = 10.0 #plasma skin depth
rpic.Te_Ti = 1.0 #T_e / T_i
rpic.c = 1.0 #Computational speed of light


rpic.oneD = False #dimensionality
rpic.twoD = True
rpic.threeD = False


#grid dimensions
rpic.Nx = 128
#rpic.Nx = 20
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

#rpic.dt = 0.10/sqrt(1.0/rpic.dx/rpic.dx)



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
#ax2.set_ylim((-15, 15))
ax2.set_ylabel(r'$E_x$')

ax3 = subplot(gs[1,1])
ax3.set_ylabel(r'$B_y$, $B_z$')

ax4 = subplot(gs[0,1])
#ax4.set_ylim((-5, 5))
ax4.set_ylabel(r'Current $J_x$')

ax5 = subplot(gs[2,0])

ax6 = fig.add_subplot(gs[2,1], projection='3d')
#ax6 = subplot(gs[2,1])



for ax in [ax1, ax2, ax4]:
    ax.set_xlim((rpic.grid_xmin, rpic.grid_xmax))
    ax.set_xlabel(r'$x$')
    ax.minorticks_on()

for ax in [ax2, ax4]:
    ax.set_xlim((-1.0, rpic.Nx))
    ax.set_xlabel(r'$x$')



#initialize particles
delgamma_i = rpic.delgamma
delgamma_e = rpic.delgamma *(rpic.mi/rpic.me)*rpic.Te_Ti



##################################################
#inject particles
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

                modulate = 0.05*np.sin(2.0*pi*x/rpic.grid_xmax)+1.0
                #modulate = 1.0
                vx, vy, vz, u = rpic.boosted_maxwellian(delgamma_e, modulate*rpic.gamma)
                #vx = 0.5
                #vy = 0.0

                #weight = q*M_macro)
                w = rpic.qe * 1.0


                cell.particles = np.concatenate((cell.particles, [[x, y, z, vx, vy, vz, w]]), axis=0) #add 6D phase space 



                #Add ion/positron to exact same place
                vxi, vyi, vzi, ui = rpic.boosted_maxwellian(delgamma_i, modulate*rpic.gamma)
                #vxi = 0.5
                #vyi = 0.0

                #wi = rpic.qe * (-1.0)
                wi = rpic.qe
                vxi *= -1.0

                #x = xmin + (xmax-xmin)*np.random.rand()
                #y = ymin + (ymax-ymin)*np.random.rand()
                #z = zmin + (zmax-zmin)*np.random.rand()

                cell.particles = np.concatenate((cell.particles, [[x, y, z, vxi, vyi, vzi, wi]]), axis=0) #add 6D phase space 



                cell.Npe += 2

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
rpic.Ex_ext[:,:,:] = 0.0


max_steps = 500
t = 0.0

#Store x-t array of a quantity
A_xt = np.zeros((rpic.Nx, max_steps))
A_xt[:,0] = rpic.Ex[:,0,0]


Nt = 50
ti = 1
tgrid = np.linspace(0.0, 24, Nt)
print tgrid

#for step in range(1, max_steps):
step = 1
while ti <= Nt:
    print " step: ",step, "/", ti-1
    print "t:", t, "/", tgrid[ti-1]

    #rpic.push_half_B()

    rpic.Vay_update_velocities(rpic.mpiGrid)

    #rpic.push_half_B()

    rpic.push_E()

    #rpic.deposit_current(rpic.mpiGrid)
    #rpic.Yee_currents()

    rpic.conserving_deposit_current(rpic.mpiGrid)
    rpic.sort_particles_between_cells(rpic.mpiGrid)


    #print "avg J=", mean(rpic.Jx)
    #print "total avg J=", sum(sum(sqrt(rpic.Jx * rpic.Jx + rpic.Jy*rpic.Jy)))

    #print "avg JY=", mean(rpic.JYx)
    #print "total avg JY=", sum(sum(sqrt(rpic.JYx * rpic.JYx + rpic.JYy*rpic.JYy)))


    #apply filters
    for sweeps in range(64):
        rpic.filter_current(0.5,1) #x sweep
        rpic.filter_current(0.5,2) #y sweep
    #rpic.filter_current(-1.0/6.0, 1) #put some power back with negative sweeps
    #rpic.filter_current(-1.0/6.0, 2)


    #I/O
    ################################################## 
    if t > tgrid[ti]:
        A_xt[:,ti] = rpic.Ex[:,0,0]
        #A_xt[:,ti] = rpic.JYx[:,0,0]


        particles = rpic.collect_grid(rpic.mpiGrid)

        electrons, positrons = rpic.divide_species(particles)

        print "shape", shape(electrons)
        x = electrons[:,0]
        y = electrons[:,1]
        z = electrons[:,2]
        vx = electrons[:,3]
        vy = electrons[:,4]
        vz = electrons[:,5]
        res1a = ax1.plot(x, vx, "k.", alpha=0.8)

        x2 = positrons[:,0]
        vx2= positrons[:,3]
        res1b = ax1.plot(x2, vx2, "b.", alpha=0.8)


        res2 = ax2.plot(XX[0,:,0], rpic.Ex[:,0,0], "b-")
        #res3 = ax3.plot(XX[0,:,0], rpic.By[:,0,0], "b-")
        #res3 = ax3.plot(XX[0,:,0], rpic.Bz[:,0,0], "r-")
        res3a = ax3.plot(rpic.JYxt[:,2,0])
        res3b = ax3.plot(rpic.JYxt[:,2,1])

        #res4a = ax4.plot(XX[0,:,0], rpic.Jx[:,0,0], "b--")
        #res4b = ax4.plot(XX[0,:,0], rpic.Jx[:,1,0], "g--")
        res4c = ax4.plot(XX[0,:,0], rpic.JYx[:,0,0], "r-")
        res4d = ax4.plot(XX[0,:,0], rpic.JYx[:,1,0], "m--")

        #ax5 = slice_pdf(ax5, electrons)

        #ax6 = plot_pdf2d(ax6, vx, vy)
        #ax6 = plot_pdf3d(ax6, vx, vy)
        #ax6 = plot_pdf3d_isosurf(ax6, vx, vy)

        fname = path+'/oneD_'+str(ti)+'.png'
        savefig(fname)

        #clean figure
        res1a.pop(0).remove()
        res1b.pop(0).remove()
        res2.pop(0).remove()
        res3a.pop(0).remove()
        res3b.pop(0).remove()
        #res4a.pop(0).remove()
        #res4b.pop(0).remove()
        res4c.pop(0).remove()
        res4d.pop(0).remove()


        ti += 1


    t += rpic.dt*rpic.wpe
    step += 1
#end of loop



