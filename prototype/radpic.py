import numpy as np
import math
from pylab import *
from mpl_toolkits.mplot3d import Axes3D



#set up figure
fig = figure(figsize=(10, 4), dpi=200)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(1, 1)
gs.update(hspace = 0.2)
gs.update(wspace = 0.2)


#physical parameters
e = 1.0
me = 1.0
c = 1.0

#pick one operating mode
oneD = True
twoD = False
threeD = False


#grid dimensions
Nx = 100
Ny = 1
Nz = 1

Nx_wrap = True
Ny_wrap = True
Nz_wrap = True


grid_xmin=0.0
grid_xmax=10.0

grid_ymin=0.0
grid_ymax=10.0

grid_zmin=0.0
grid_zmax=10.0

Np = 1000 #total number of particles


#derived values
xgrid = np.linspace(grid_xmin, grid_xmax, Nx+1)
ygrid = np.linspace(grid_ymin, grid_ymax, Ny+1)
zgrid = np.linspace(grid_zmin, grid_zmax, Nz+1)

#Mesh grid for plotting
XX, YY, ZZ = np.meshgrid(np.linspace(grid_xmin, grid_xmax, Nx),
                         np.linspace(grid_ymin, grid_ymax, Ny),
                         np.linspace(grid_zmin, grid_zmax, Nz)
                        )


Ncells = Nx*Ny*Nz

dx = diff(xgrid)[0]
dy = diff(ygrid)[0]
dz = diff(zgrid)[0]


dt = 0.99/sqrt(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz)

#correct grid sizes
if oneD:
    Ny = 1
    Nz = 1
if twoD:
    Nz = 1


#initialize B, E, and J fields (in global scope)
Bx = np.zeros((Nx,Ny,Nz))
By = np.zeros((Nx,Ny,Nz))
Bz = np.zeros((Nx,Ny,Nz))

Ex = np.zeros((Nx,Ny,Nz))
Ey = np.zeros((Nx,Ny,Nz))
Ez = np.zeros((Nx,Ny,Nz))

Jx = np.zeros((Nx,Ny,Nz))
Jy = np.zeros((Nx,Ny,Nz))
Jz = np.zeros((Nx,Ny,Nz))



class CellClass(object):
    def __init__(self):
        self.Npe = 0
        self.electrons = np.empty((0,6), dtype=float64)        

        #self.xmin = 0
        #self.xmax = 0
        #self.ymin = 0
        #self.ymax = 0
        #self.zmin = 0
        #self.zmax = 0

def grid_limits(i,j,k):
    xmin = xgrid[i]
    xmax = xgrid[i+1]

    ymin = ygrid[j]
    ymax = ygrid[j+1]

    zmin = zgrid[k]
    zmax = zgrid[k+1]

    return xmin,xmax,ymin,ymax,zmin,zmax

def grid_lengths(i,j,k):
    dx = xgrid[i+1]-xgrid[i]
    dy = ygrid[j+1]-ygrid[j]
    dz = zgrid[k+1]-zgrid[k]

    return dx,dy,dz


#draw samples from Maxwellian distribution using rejection sampling
def Maxwellian(vb):
    fmax = 0.5*(1.0 + exp(-2.0*vb*vb))
    vmin = -5.0*vb
    vmax =  5.0*vb
    vf = vmin + (vmax-vmin)*np.random.rand()

    f = 0.5*(exp(-(vf-vb)*(vf-vb)/2.0))

    x = fmax*np.random.rand()

    if x > f:
        return Maxwellian(vb)

    return vf




#create grid
mpiGrid = np.empty((Nx,Ny,Nz), dtype=np.object)



#inject particles
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            cell = CellClass() 

            #cell limits
            xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
            
            #Temperature
            vb = 2.0

            for n in range(Np/Ncells):
                #random uniform location inside cell
                x = xmin + (xmax-xmin)*np.random.rand()
                y = ymin + (ymax-ymin)*np.random.rand()
                z = zmin + (zmax-zmin)*np.random.rand()

                vx = Maxwellian(vb)
                vy = Maxwellian(vb)
                vz = Maxwellian(vb)

                cell.electrons = np.concatenate((cell.electrons, [[x, y, z, vx, vy, vz]]), axis=0) #add 6D phase space 

                cell.Npe += 1

            mpiGrid[i,j,k] = cell



#create E field
#Ex[:,:,:] = 0.1
#Ey[:,:,:] = 0.1
#Ez[:,:,:] = 0.1



#deposit particle currents into the mesh
def deposit_current(grid):
    global e
    global me
    global me

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                cell = grid[i,j,k]
                particles = cell.electrons

                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                rhop = e/(dx*dy*dz)
                q = 1.0

                x = particles[:,0]
                y = particles[:,1]
                z = particles[:,2]
                ux = particles[:,3]
                uy = particles[:,4]
                uz = particles[:,5]

                gamma = sqrt(1.0 + ux*ux + uy*uy + uz*uz)
                ux = ux*c/gamma
                uy = uy*c/gamma
                uz = uz*c/gamma

                fp = (x - xmin)/dx
                fq = (y - ymin)/dy
                fr = (z - zmin)/dz

                dJx = 0.5*sum(ux)*q*rhop
                dJy = 0.5*sum(uy)*q*rhop
                dJz = 0.5*sum(uz)*q*rhop
            
                #TODO implement higher order cloud-in-the-cell models
                Jx[i,j,k] = dJx
                Jy[i,j,k] = dJy
                Jz[i,j,k] = dJz

    return


#Vay pusher and particle propagator
# note: we update every particle at once in vector format
def update_velocities(grid):
    global e
    global me
    global dt

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                cell = grid[i,j,k]
                particles = cell.electrons

                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                q = 1.0
                m = me

                x = particles[:,0]
                y = particles[:,1]
                z = particles[:,2]
                ux = particles[:,3]
                uy = particles[:,4]
                uz = particles[:,5]

                xd = (x - xmin)/dx
                yd = (y - xmin)/dy
                zd = (z - xmin)/dz

                #interpolate E and B
                # now we just pick the grid value directly
                Exi = Ex[i,j,k]
                Eyi = Ey[i,j,k]
                Ezi = Ez[i,j,k]
                
                Bxi = Bx[i,j,k]
                Byi = By[i,j,k] 
                Bzi = Bz[i,j,k] 

                uxm = ux + q*e*Exi*dt/(2.0*m*c)
                uym = uy + q*e*Eyi*dt/(2.0*m*c)
                uzm = uz + q*e*Ezi*dt/(2.0*m*c)


                #Lorentz transform
                gamma = sqrt(1.0 + uxm*uxm + uym*uym + uzm*uzm)
                #gamma = 1.0

                #Calculate u'
                tx = q*e*Bxi*dt/(2.0*gamma*m*c)
                ty = q*e*Byi*dt/(2.0*gamma*m*c)
                tz = q*e*Bzi*dt/(2.0*gamma*m*c)

                ux0 = uxm + uym*tz - uzm*ty
                uy0 = uym + uzm*tx - uxm*tz
                uz0 = uzm + uxm*ty - uym*tx

                #calculate u+
                sx = 2.0*tx/(1.0 + tx*tx + ty*ty + tz*tz)
                sy = 2.0*ty/(1.0 + tx*tx + ty*ty + tz*tz)
                sz = 2.0*tz/(1.0 + tx*tx + ty*ty + tz*tz)

                uxp = uxm + uy0*sz - uz0*sy
                uyp = uym + uz0*sx - ux0*sz
                uzp = uzm + ux0*sy - uy0*sx

                # update velocities
                #t(dt/2) -> t(dt/2)
                particles[:,3] = uxp + q*e*Exi*dt/(2.0*m*c)
                particles[:,4] = uyp + q*e*Eyi*dt/(2.0*m*c)
                particles[:,5] = uzp + q*e*Ezi*dt/(2.0*m*c)
                
                #update locations and propagate particles
                uxn = particles[:,3]
                uyn = particles[:,4]
                uzn = particles[:,5]
                
                gamma = sqrt(1.0 + uxn*uxn + uyn*uyn + uzn*uzn)
                #gamma = 1.0

                if oneD:
                    particles[:,0] = particles[:,0] + (c*dt/gamma)*uxn
                if twoD:
                    particles[:,0] = particles[:,0] + (c*dt/gamma)*uxn
                    particles[:,1] = particles[:,1] + (c*dt/gamma)*uyn
                if threeD:
                    particles[:,0] = particles[:,0] + (c*dt/gamma)*uxn
                    particles[:,1] = particles[:,1] + (c*dt/gamma)*uyn
                    particles[:,2] = particles[:,2] + (c*dt/gamma)*uzn


    return


# sort particles between neighboring grid cells
def sort_particles_between_cells(grid):

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):

                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                #print "cell=",i,j,k, "x=",xmin,xmax,dx
                #print "cell=",i,j,k, "y=",ymin,ymax,dy
                #print "cell=",i,j,k, "z=",zmin,zmax,dz
                #print "x=", x[0], "new x:", round((x[0]-xmin)/dx - 1)
                #print "y=", y[0], "new y:", round((y[0]-ymin)/dy - 1)
                #print "z=", z[0], "new z:", round((z[0]-zmin)/dz - 1)

                ip = 0
                while (ip < grid[i,j,k].Npe):
                    particle = grid[i,j,k].electrons[ip,:]
                    iddf = (particle[0] - xmin - dx)/dx
                    jddf = (particle[1] - ymin - dy)/dy
                    kddf = (particle[2] - zmin - dz)/dz

                    idd = int(iddf) if iddf < 0 else int(iddf+1)
                    jdd = int(jddf) if jddf < 0 else int(jddf+1)
                    kdd = int(kddf) if kddf < 0 else int(kddf+1)

                    #print "particle coords=", grid[i,j,k].electrons[ip,:]
                    #print "ip=",ip, "indx=", idd, jdd, kdd

                    if idd == 0 and jdd == 0 and kdd == 0:
                        #print "inside"
                        ip += 1
                    else:
                        newi = i + idd
                        newj = j + jdd
                        newk = k + kdd

                        #periodic boundary conditions
                        remove_particle = False
                        #print "before wrap=",newi, newj, newk

                        if newi < 0:
                            if Nx_wrap:
                                newi += Nx
                                particle[0] += grid_xmax
                            else:
                                remove_particle = True
                        if newj < 0:
                            if Ny_wrap:
                                newj += Ny
                                particle[1] += grid_ymax
                            else:
                                remove_particle = True
                        if newk < 0:
                            if Nz_wrap:
                                newk += Nz
                                particle[2] += grid_zmax
                            else:
                                remove_particle = True

                        if newi >= Nx:
                            if Nx_wrap:
                                newi -= Nx
                                particle[0] -= grid_xmax
                            else:
                                remove_particle = True
                        if newj >= Ny:
                            if Ny_wrap:
                                newj -= Ny
                                particle[1] -= grid_ymax
                            else:
                                remove_particle = True
                        if newk >= Nz:
                            if Nz_wrap:
                                newk -= Nz
                                particle[2] -= grid_zmax
                            else:
                                remove_particle = True


                        #add to new cell
                        if (not remove_particle):
                            #print "shifting particle", newi, newj, newk
                            #print "new population"
                            #print grid[newi,newj,newk].electrons

                            grid[newi, newj, newk].electrons = np.concatenate( (grid[newi,newj,newk].electrons, [particle] ), axis=0)
                            grid[newi, newj, newk].Npe += 1

                            #print "after populating"
                            #print grid[newi,newj,newk].electrons

                        #remove from previous cell
                        np.delete(grid[i,j,k].electrons, ip, 0)
                        grid[i,j,k].Npe -= 1
                        ip += 1

    return 


def push_half_B():
    EY = np.zeros(3)
    EYp1 = np.zeros(3)

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                EY[0] = Ex[i,j,k]                
                EY[1] = Ey[i,j,k]                
                EY[2] = Ez[i,j,k]                
                
                ii = i+1 if (i < Nx 
                jj = j+1
                kk = k+1
                EYp1[0] = Ex[ii,jj,kk]                
                EYp1[1] = Ey[ii,jj,kk]                
                EYp1[2] = Ez[ii,jj,kk]                



    return


##################################################

def collect_grid(grid):
    electrons = np.empty((0,6), dtype=float64)

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                cell = grid[i,j,k]
                electrons = np.concatenate((electrons, cell.electrons), axis=0)
    return electrons


#####


deposit_current(mpiGrid)


if oneD:
    ax1 = subplot(gs[0,0])
    ax1.set_xlim((grid_xmin,grid_xmax))
    ax1.set_ylim((-5, 5))
if twoD:
    ax1 = subplot(gs[0,0])
    ax1.set_xlim((grid_xmin,grid_xmax))
    ax1.set_ylim((grid_ymin,grid_ymax))
if threeD:
    fig = plt.figure()
    ax1 = fig.add_subplot(224, projection='3d')

    ax2 = fig.add_subplot(221, projection='3d')
    ax3 = fig.add_subplot(222, projection='3d')
    ax4 = fig.add_subplot(223, projection='3d')

    ax1.set_xlim((grid_xmin,grid_xmax))
    ax1.set_ylim((grid_ymin,grid_ymax))
    ax1.set_zlim((grid_zmin,grid_zmax))

    ax2.set_title(r'$\bar{E}$')
    ax2.set_xlim((grid_xmin,grid_xmax))
    ax2.set_ylim((grid_ymin,grid_ymax))
    ax2.set_zlim((grid_zmin,grid_zmax))

    ax3.set_title(r'$\bar{B}$')
    ax3.set_xlim((grid_xmin,grid_xmax))
    ax3.set_ylim((grid_ymin,grid_ymax))
    ax3.set_zlim((grid_zmin,grid_zmax))

    ax4.set_title(r'Current $\bar{J}$')
    ax4.set_xlim((grid_xmin,grid_xmax))
    ax4.set_ylim((grid_ymin,grid_ymax))
    ax4.set_zlim((grid_zmin,grid_zmax))



##################################################
##################################################
##################################################


max_steps = 5

for step in range(1, max_steps):
    print " step: ",step

    #push_half_B

    update_velocities(mpiGrid)
    
    sort_particles_between_cells(mpiGrid)

    #push half B
    #push E

    deposit_current(mpiGrid)
    #Yee_currents(mpiGrid)

    #apply filters

    #I/O
    electrons = collect_grid(mpiGrid)
    x = electrons[:,0]
    y = electrons[:,1]
    z = electrons[:,2]
    vx = electrons[:,3]
    vy = electrons[:,4]
    vz = electrons[:,5]

    if oneD:
        res1 = ax1.plot(x, vx, "k.")

        savefig('pic1d_'+str(step)+'.png')
        res1.pop(0).remove()
    if twoD:
        res1 = ax1.plot(x, y, "k.")


        savefig('pic2d_'+str(step)+'.png')
        res1.pop(0).remove()

    if threeD:
        res1  = ax1.plot(x, y, z, "k.")

        res2 =  ax2.quiver(XX, YY, ZZ, Ex, Ey, Ez, pivot='tail')

        res3 =  ax3.quiver(XX, YY, ZZ, Bx, By, Bz, pivot='tail')

        res4 =  ax4.quiver(XX, YY, ZZ, Jx, Jy, Jz, pivot='tail')

        savefig('pic3d_'+str(step)+'.png')

        #clear figure
        res1.pop(0).remove()
        ax2.collections.remove(res2)
        ax3.collections.remove(res3)
        ax4.collections.remove(res4)


#end of loop








