import numpy as np
import math
from pylab import *

#set seed to get reproducible errors & results
np.random.seed(0)

#Computational units
Nppc = 0
gamma = 0.0
delgamma = 0.0
me = 1.0
mi = 1.0
qe = 1.0
qi = 1.0
qme = 1.0
qmi = 1.0
delta = 0.0
wpe = 0.0
c = 1.0


#plasmaFreq = 0.0
#skinDepth = 0.0


#pick one operating mode
zeroD = False
oneD = False
twoD = False
threeD = False


#grid dimensions
Nx = 1
Ny = 1
Nz = 1

Nx_wrap = False
Ny_wrap = False
Nz_wrap = False


grid_xmin=0.0
grid_xmax=1.0

grid_ymin=0.0
grid_ymax=1.0

grid_zmin=0.0
grid_zmax=1.0
Np = 0 #total number of particles


#derived values
xgrid = np.linspace(grid_xmin, grid_xmax, Nx+1)
ygrid = np.linspace(grid_ymin, grid_ymax, Ny+1)
zgrid = np.linspace(grid_zmin, grid_zmax, Nz+1)

Ncells = Nx*Ny*Nz

dx = diff(xgrid)[0]
dy = diff(ygrid)[0]
dz = diff(zgrid)[0]

dt = 0.0

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

JYxt = np.zeros((Nx+2,Ny+2,Nz+2))
JYyt = np.zeros((Nx+2,Ny+2,Nz+2))
JYzt = np.zeros((Nx+2,Ny+2,Nz+2))

JYx = np.zeros((Nx,Ny,Nz))
JYy = np.zeros((Nx,Ny,Nz))
JYz = np.zeros((Nx,Ny,Nz))


#external fields
Bx_ext = np.zeros((Nx,Ny,Nz))
By_ext = np.zeros((Nx,Ny,Nz))
Bz_ext = np.zeros((Nx,Ny,Nz))

Ex_ext = np.zeros((Nx,Ny,Nz))
Ey_ext = np.zeros((Nx,Ny,Nz))
Ez_ext = np.zeros((Nx,Ny,Nz))

#radiation fields

#From Zeltron
Uph=0.0 #external photon energy density

#synchrotron energy losses
Esyned=0.0 
Esynpd=0.0 
Esyneb=0.0 
Esynpb=0.0 
Esyn_electrons=0.0 
Esyn_ions=0.0 

#inverse Compton energy losses
Eicsed=0.0 
Eicspd=0.0 
Eicseb=0.0 
Eicspb=0.0 
Eics_electrons=0.0 
Eics_ions=0.0 

#external radiation field energy density
B0 = 0.0
udens_ratio = 1.0
Uph = udens_ratio*(B0*B0/(8.0*pi))

#em_energy
#synchrotrhon
#rad_energy(esyneb, esyned, esyn_electrons, electrons, syn)
#rad_energy(esynpb, esynpd, esyn_ions, ions, syn)

#inverse compton
#rad_energy(eicseb, eicsed, eics_electrons, electrons, syn)
#rad_energy(eicspb, eicspd, eics_ions, ions, syn)

#spectrum and angular distribution
#this subroutine computes particles spectrum and angular distribution
#input 
#pcl #particel distribution function
#it timestep
#spec particle spectrum
#name of particle
#output:
# u, phi, lmabda, dN/dOmega/du, dNdu


#spectrum_angular(pcl_ed, electrons, drift #drift electrons
#spectrum_angular(pcl_pd, ions, drift #drift ions

#spectrum_angular(pcl_ed, electrons, bg #bkg electrons
#spectrum_angular(pcl_pd, ions, bg #bkg ions

#analyze radiation
# Computes total sync. radiatin spectrum and angular disteibution
# input: 
#  pcl particle distribution
#  Bxg x-component of B at nodes at t
#  Byg y-component of B at nodes at t
#  Bzg z-component of B at nodes at t
#  it timestep
#  spec particle species
#  sym name of the particles
# output nu(Hz), nuFnu(erg/s)o
#analysis_sync(mi, pclpb, Bxg, Byg, Bzg...


#Zeltron Boris pusher
# #Esyn total sync. energy losses between t and t+dt
# #Eics total inverse Compton energy losses between t and t+dt
# #external photon energy density



#create grid
mpiGrid = np.empty((Nx,Ny,Nz), dtype=np.object)

#get lengths from skin depths
def grid_scale_skindepth(delta):
    global grid_xmin, grid_xmax, Nx
    global grid_ymin, grid_ymax, Ny
    global grid_zmin, grid_zmax, Nx
    global dx, dy, dz
    
    #start from origo
    grid_xmin = 0.0
    grid_ymin = 0.0
    grid_zmin = 0.0

    dx = (delta**2) *(e*e/(c*c*me)) * Nppc
    
    grid_xmax = Nx*dx
    grid_ymax = Ny*dx
    grid_zmax = Nz*dx

    # try setting d elements instead
    dy = dx
    dz = dx

    return



def init():

    #fix dimensionality
    global zeroD
    global oneD
    global twoD
    global threeD
    global Nx
    global Ny
    global Nz

    if zeroD:
        oneD = False
        twoD = False
        threeD = False
    elif oneD:
        twoD = False
        threeD = False
    elif twoD:
        oneD = False
        threeD = False
    elif threeD:
        oneD = False
        twoD = False


    #correct grid sizes
    if zeroD:
        Nx = 1
        Ny = 1
        Nz = 1
    if oneD:
        Ny = 1
        Nz = 1
    if twoD:
        Nz = 1

    global Ncells
    Ncells = Nx*Ny*Nz



    #plasma frequency
    #global wpe, c, delta
    #wpe = c/delta
    #print "Plasma frequency w_pe:", wpe
    
    #if gamma < 1 it is interpreted as v/c
    #global gamma
    #if gamma < 1.0:
    #    print "non rel. v/c:", gamma 
    #    gamma = sqrt(1.0/(1.0 - gamma**2))
    #print "rel. gamma:", gamma

    #drift velocity
    #global beta
    #beta = sqrt(1.0 - 1.0/gamma**2)
    #print "Drift velocity:", beta

    #electron charge in computational units
    #global qe, me, mi, Nppc
    #qe = -(gamma*wpe**2)/((Nppc*0.5) * (1.0 + me/mi))
    #print "electron charge in computational units:", qe

    #ion/positron charge in computational units
    #global qi
    #qi = - qe
    #print "ion charge in computational units:", qi

    #electron mass
    #me = me*abs(qi)
    #print "electron mass:", me

    #ion/positron mass
    #mi = mi*abs(qi)
    #print "ion mass:", mi
    

    #charge to mass ratios
    #global qme, qmi
    #qme = qe/me #equals -1/me in input
    #qmi = qi/mi #equals -1/me in input
    #print "electron/ion charge to mass ratios:", qme, qmi
    
    
    #Binit
    #Binit = sqrt(gamma*ppc0*0.5*c**2*(me*(1.0+me/mi))*sigma)
    

    #computational units
    #global dx
    #global dy
    #global dz
    #dx = 1.0
    #dy = 1.0
    #dz = 1.0

    #global xgrid
    #global ygrid
    #global zgrid
    #xgrid = np.linspace(grid_xmin, grid_xmax, Nx+1)
    #ygrid = np.linspace(grid_ymin, grid_ymax, Ny+1)
    #zgrid = np.linspace(grid_zmin, grid_zmax, Nz+1)


    global xgrid
    global ygrid
    global zgrid
    
    global dx
    global dy
    global dz

    xgrid = np.linspace(grid_xmin, grid_xmax, Nx+1)
    ygrid = np.linspace(grid_ymin, grid_ymax, Ny+1)
    zgrid = np.linspace(grid_zmin, grid_zmax, Nz+1)
    #xgrid = np.arange(grid_xmin, grid_xmax, dx)
    #ygrid = np.arange(grid_ymin, grid_ymax, dy)
    #zgrid = np.arange(grid_zmin, grid_zmax + dz, dz)

    dx = diff(xgrid)[0]
    dy = diff(ygrid)[0]
    dz = diff(zgrid)[0]

    global dt 
    dt = 0.45/sqrt(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz)


    #initialize B, E, and J fields (in global scope)
    global Bx
    global By
    global Bz
    global Ex
    global Ey
    global Ez
    global Jx
    global Jy
    global Jz
    global JYx
    global JYy
    global JYz

    global Bx_ext
    global By_ext
    global Bz_ext
    global Ex_ext
    global Ey_ext
    global Ez_ext

    Bx = np.zeros((Nx,  Ny, Nz))
    By = np.zeros((Nx,  Ny, Nz))
    Bz = np.zeros((Nx,  Ny, Nz))
    
    Bx_ext = np.zeros((Nx,  Ny, Nz))
    By_ext = np.zeros((Nx,  Ny, Nz))
    Bz_ext = np.zeros((Nx,  Ny, Nz))

    Ex = np.zeros((Nx,  Ny, Nz))
    Ey = np.zeros((Nx,  Ny, Nz))
    Ez = np.zeros((Nx,  Ny, Nz))

    Ex_ext = np.zeros((Nx,  Ny, Nz))
    Ey_ext = np.zeros((Nx,  Ny, Nz))
    Ez_ext = np.zeros((Nx,  Ny, Nz))
    
    Jx = np.zeros((Nx,  Ny, Nz))
    Jy = np.zeros((Nx,  Ny, Nz))
    Jz = np.zeros((Nx,  Ny, Nz))
    
    JYx = np.zeros((Nx, Ny, Nz))
    JYy = np.zeros((Nx, Ny, Nz))
    JYz = np.zeros((Nx, Ny, Nz))
    
    global JYxt, JYyt, JYzt
    JYxt = np.zeros((Nx+2,Ny+2,Nz+2))
    JYyt = np.zeros((Nx+2,Ny+2,Nz+2))
    JYzt = np.zeros((Nx+2,Ny+2,Nz+2))

    #create grid
    global mpiGrid
    mpiGrid = np.empty((Nx, Ny, Nz), dtype=np.object)

    #print "##################################################"
    print "Nx Ny Nz", Nx, Ny, Nz 
    print "xmin xmax", grid_xmin, grid_xmax
    print "ymin ymax", grid_ymin, grid_ymax
    print "zmin zmax", grid_zmin, grid_zmax
    #print "dt", dt


    global ne, wpe, kT, r_deb, v_eth

    #number density
    ne_x = Np/(grid_xmax - grid_xmin)
    ne_xy = Np/(grid_xmax - grid_xmin)/(grid_ymax - grid_ymin)
    ne   = Np/(grid_xmax - grid_xmin)/(grid_ymax - grid_ymin)/(grid_zmax - grid_zmin)

    ##plasma frequency
    wpe = sqrt((ne*e**2)/me)
    #wpe = sqrt((ne*e**2)/me)

    ##plasma skin depth
    delta = c/wpe

    ##Debye length
    r_deb = sqrt(delgamma/ne_xy)

    ##electron thermal velocity
    v_eth = sqrt(delgamma/me)

    print "ne (1d):", ne_x
    print "ne (2d):", ne_xy
    print "ne (3d):", ne
    print "wpe:", wpe, "tpw=", 1.0/wpe
    print "dt/tpw =", dt*wpe

    print "delta x:", delta, delta/dx
    print "delta y:", delta, delta/dy
    print "delta z:", delta, delta/dz

    print "r_deb x:", r_deb, r_deb/dx
    print "r_deb y:", r_deb, r_deb/dy

    print "v_eth:", v_eth


    return

##################################################


class CellClass(object):
    def __init__(self):
        self.Npe = 0
        self.particles = np.empty((0,7), dtype=float64)        



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

def xbc(i):
    #return i+1
    if i < 0:
        return int(i + Nx)
    if i >= Nx:
        return int(i - Nx)
    return int(i)

def ybc(j):
    #return j+1
    if j < 0:
        return int(j + Ny)
    if j >= Ny:
        return int(j - Ny)
    return int(j)

def zbc(k):
    if k < 0:
        return int(k + Nz)
    if k >= Nz:
        return int(k - Nz)
    return int(k)


#draw samples from Maxwellian distribution using rejection sampling
def Maxwellian(vb, theta):
    fmax = 0.5*(1.0 + exp(-vb*vb/(2.0*theta*theta)))
    vmin = -5.0*vb
    vmax =  5.0*vb
    vf = vmin + (vmax-vmin)*np.random.rand()

    f = 0.5*(exp(-(vf-vb)*(vf-vb) / (2.0*theta*theta)))

    x = fmax*np.random.rand()

    if x > f:
        return Maxwellian(vb, theta)

    return vf



#deposit particle currents into the mesh
def deposit_current(grid):
    global qe
    global me

    np1 = 0
    np2 = 0
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                cell = grid[i,j,k]
                particles = cell.particles
                
                #print "cell:",i,j,k, " N=", len(particles), "/", cell.Npe
                np1 += len(particles)
                np2 += cell.Npe


                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                x = particles[:,0]
                y = particles[:,1]
                z = particles[:,2]
                ux = particles[:,3]
                uy = particles[:,4]
                uz = particles[:,5]
            
                #rhop = qe/(dx*dy*dz)
                #q = 1.0
                #full particle weight vector, i.e. we push all of the particles at the same time
                #rhop = particles[:,6]/(dx*dy*dz) #w contains q*m
                rhop = particles[:,6]/(dx*dy) #w contains q*m
                #rhop = particles[:,6] #charge density assuming dx=dy=dz = 1

                gamma = sqrt(1.0 + ux*ux + uy*uy + uz*uz)
                ux = ux*c/gamma
                uy = uy*c/gamma
                uz = uz*c/gamma

                fp = (x - xmin)/dx
                fq = (y - ymin)/dy
                fr = (z - zmin)/dz

            
                #just get current to center cell; 0th order model
                if zeroD:
                    dJx = 0.5*sum(ux*rhop)
                    dJy = 0.5*sum(uy*rhop)
                    dJz = 0.5*sum(uz*rhop)
                    Jx[i,j,k] += dJx
                    Jy[i,j,k] += dJy
                    Jz[i,j,k] += dJz

                #cloud-in-the-cell model
                if oneD:
                    for xdir in [0,1]:
                        wx = 1.0-fp if xdir == 0 else fp
                        ii = i+xdir
                        if ii >= Nx:
                            ii -= Nx
                        dJx = 0.5*sum(ux*wx*rhop)
                        dJy = 0.5*sum(uy*wx*rhop)
                        dJz = 0.5*sum(uz*wx*rhop)
                        Jx[ii,j,k] += dJx
                        Jy[ii,j,k] += dJy
                        Jz[ii,j,k] += dJz

                if twoD:
                    for xdir in [0,1]:
                        wx = 1.0-fp if xdir == 0 else fp
                        ii = i+xdir
                        if ii >= Nx:
                            ii -= Nx
                        for ydir in [0,1]:
                            wy = 1.0-fq if ydir == 0 else fq
                            jj = j+ydir
                            if jj >= Ny:
                                jj -= Ny

                            dJx = 0.5*sum(ux*wx*wy*rhop)
                            dJy = 0.5*sum(uy*wx*wy*rhop)
                            dJz = 0.5*sum(uz*wx*wy*rhop)
                            Jx[ii,jj,k] += dJx
                            Jy[ii,jj,k] += dJy
                            Jz[ii,jj,k] += dJz

                if threeD:
                    for xdir in [0,1]:
                        wx = 1.0-fp if xdir == 0 else fp
                        ii = i+xdir
                        if ii >= Nx:
                            ii -= Nx
                        for ydir in [0,1]:
                            wy = 1.0-fq if ydir == 0 else fq
                            jj = j+ydir
                            if jj >= Ny:
                                jj -= Ny
                            for zdir in [0,1]:
                                wz = 1.0-fr if zdir == 0 else fr
                                kk = k+zdir
                                if kk >= Nz:
                                    kk -= Nz

                                dJx = 0.5*sum(ux*wx*wy*wz*rhop)
                                dJy = 0.5*sum(uy*wx*wy*wz*rhop)
                                dJz = 0.5*sum(uz*wx*wy*wz*rhop)
                                Jx[ii,jj,kk] += dJx
                                Jy[ii,jj,kk] += dJy
                                Jz[ii,jj,kk] += dJz



    #print "    TOTAL=", np1, "/", np2
    return

#deposit particle currents into the mesh
def conserving_deposit_current(grid):
    global qe
    global me

    global JYx, JYy, JYz

    global JYxt, JYyt, JYzt
    
    #clean currents
    JYx[:,:,:] = 0.0
    JYy[:,:,:] = 0.0
    JYz[:,:,:] = 0.0

    JYxt[:,:,:] = 0.0
    JYyt[:,:,:] = 0.0
    JYzt[:,:,:] = 0.0

    np1 = 0
    np2 = 0
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                cell = grid[i,j,k]
                particles = cell.particles
                
                #print "cell:",i,j,k, " N=", len(particles), "/", cell.Npe
                np1 += len(particles)
                np2 += cell.Npe


                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                x2 = particles[:,0]
                y2 = particles[:,1]
                z2 = particles[:,2]
                ux = particles[:,3]
                uy = particles[:,4]
                uz = particles[:,5]
            
                #charge density vector
                #i.e. we push all of the particles at the same time
                #rhop = particles[:,6]/(dx*dy*dz) #w contains q*m

                #go back in time -dt and see what was the starting position
                #FIXME could we not just save this?
                invgamma = 1.0/sqrt(1.0 + ux*ux + uy*uy + uz*uz)
                x1 = x2 - ux*invgamma*c*dt
                y1 = y2 - uy*invgamma*c*dt
                z1 = z2 - uz*invgamma*c*dt

                q = particles[:,6]*qe 

                #x1sp = x0
                #x2sp = x
                #y1sp = y0
                #y2sp = y
                #z1=z0 #FIXME what is this?
                #z2=z

                #fp = (x - xmin)/dx
                #fq = (y - ymin)/dy
                #fr = (z - zmin)/dz

                #i1 = np.floor(x1sp/dx)
                i1 = np.floor((x1 - xmin)/dx) + i
                #i1 = np.ones(cell.Npe)*i
                #i2 = np.floor((x2)/dx)
                i2 = np.ones(cell.Npe)*i
                #i2 = i
                #j1 = np.floor(y1sp/dy)
                j1 = np.floor((y1 - ymin)/dy) + j
                #j1 = np.ones(cell.Npe)*j
                #j2 = np.floor((y2)/dy)
                #j2tmp = np.floor(y2/dy)
                j2 = np.ones(cell.Npe)*j
                #j2 = j
                #k1 = int(z1)
                #k2 = int(z2)

                xrr = np.empty_like(x1)
                yrr = np.empty_like(y1)

                for ppp in range(cell.Npe):
                    #print "---:", ppp,"  / ",i,j,k
                    #print "x1=", x1[ppp], " y1=", y1[ppp], " z1=",z1[ppp]
                    #print "x2=", x2[ppp], " y2=", y2[ppp], " z2=",z2[ppp]
                    #print "i1 / i2", i1[ppp], i2[ppp]
                    #print "j1 / j2", j1[ppp], j2[ppp]

                    #if i1[ppp] != i:
                    #    print "ERRRRRRRR X", i1[ppp], i
                    #if j1[ppp] != j:
                    #    print "ERRRRRRRR Y", j1[ppp], j

                    #xr = 0.0
                    #if i1[ppp] == i2:
                    #    print "i1 = i2;",
                    #    xr = 0.5*(x1[ppp] + x2[ppp])
                    #else:
                    #    print " != ", max(i1[ppp], i2)
                    #    xri = max(i1[ppp], i2)
                    #    xr = xgrid[xri]
                    #print "xr=", xr

                    i1p = int(i1[ppp])
                    i2p = int(i2[ppp])
                    xr2 = min(
                            min(xgrid[i1p], xgrid[i2p]) + dx,
                            max(max(xgrid[i1p], xgrid[i2p]), 0.5*(x1[ppp] + x2[ppp]))
                             )
                    #print "xr minmax=", xr2
                    xrr[ppp] = xr2

                    j1p = int(j1[ppp])
                    j2p = int(j2[ppp])
                    yr2 = min(
                            min(ygrid[j1p], ygrid[j2p]) + dy,
                            max(max(ygrid[j1p], ygrid[j2p]), 0.5*(y1[ppp] + y2[ppp]))
                             )
                    #print "yr minmax=", yr2
                    yrr[ppp] = yr2
                #print "end of cell"
            
                #build currents; 
                # not including -q to be consistent with notation
                Fx1 = -q*(xrr - x1)/dt
                Fy1 = -q*(yrr - y1)/dt

                Fx2 = -q*(x2 - xrr)/dt
                Fy2 = -q*(y2 - yrr)/dt

                Wx1 = 0.5*(x1 + xrr)/dx -i1
                Wy1 = 0.5*(y1 + yrr)/dy -j1

                Wx2 = 0.5*(xrr + x2)/dx -i2
                Wy2 = 0.5*(yrr + y2)/dy -j2

                #print "Fx1", Fx1
                #print "Fy1", Fy1
                #print "Fx2", Fx2
                #print "Fy2", Fy2
                #print "Wx1", Wx1
                #print "Wy1", Wy1
                #print "Wx2", Wx2
                #print "Wy2", Wy2

                #x-dir
                Jx1a = Fx1*(1.0 - Wy1)/dx/dy
                Jx1b = Fx1*Wy1/dx/dy

                Jx2a = Fx2*(1.0 - Wy2)/dx/dy
                Jx2b = Fx2*Wy2/dx/dy
                
                #y-dir
                Jy1a = Fy1*(1.0 - Wx1)/dx/dy
                Jy1b = Fy1*Wx1/dx/dy

                Jy2a = Fy2*(1.0 - Wx2)/dx/dy
                Jy2b = Fy2*Wx2/dx/dy
                
                for ppp in range(cell.Npe):
                    #printflag=False
                    #if xbc(i1[ppp]) == 5:
                    #    printflag=True
                    #    print "--- i1", i1[ppp], "j1", j1[ppp] 

                    i1p = xbc(i1[ppp])
                    j1p = ybc(j1[ppp])

                    i2p = xbc(i2[ppp])
                    j2p = ybc(j2[ppp])
                    
                    i1p1 = xbc(i1p + 1)
                    j1p1 = ybc(j1p + 1)

                    i2p1 = xbc(i2p + 1)
                    j2p1 = ybc(j2p + 1)

                    #if i1p == 0 or i1p1 == 0 or i2p1 == 0 or i2 == 0:

                    #if i1p == Nx-1:
                    #    print "i1p 0"
                    #if i1p1 == Nx-1:
                    #    print "i1p1 0"
                    #if i2 == Nx-1:
                    #    print "i2 0"
                    #print "i",i1p, " j", j1p1, "k",k

                    #if printflag:
                    #    print "   i1p", i1p, "j1p", j1p, "J", Jx1a[ppp]

                    JYx[i1p, j1p, k] += Jx1a[ppp]
                    #JYx[i1p, j1p, k] += 1.0
                    JYx[i1p, j1p1,k] += Jx1b[ppp]
                    JYx[i2p, j2p,  k] += Jx2a[ppp]
                    JYx[i2p, j2p1,k] += Jx2b[ppp]

                    JYy[i1p, j1p, k] += Jy1a[ppp] 
                    JYy[i1p1,j1p, k] += Jy1b[ppp]
                    JYy[i2p, j2p, k] += Jy2a[ppp]
                    JYy[i2p1,j2p, k] += Jy2b[ppp]

    #for k in range(Nz):
    #    for j in range(Ny):
    #        for i in range(Nx):
    #            JYx[i,j,k] = JYxt[i+1, j+1, k+1]
    #            JYy[i,j,k] = JYyt[i+1, j+1, k+1]
             

    return


#Yee shift current vector into staggered grid
def Yee_currents():

    Jp1 = np.zeros(3)
    J0 = np.zeros(3)

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                J0[0] = Jx[i,j,k]
                J0[1] = Jy[i,j,k]
                J0[2] = Jz[i,j,k]

                #Boundary  conditions and wrapping
                # TODO wrapping in every dimension done automatically here
                # correct this to check for Nx/y/z_wrap

                if oneD:
                    ii = i+1 if i < Nx-1 else 0
                    jj = j
                    kk = k
                elif twoD:
                    ii = i+1 if i < Nx-1 else 0
                    jj = j+1 if j < Ny-1 else 0
                    kk = k
                elif threeD:
                    ii = i+1 if i < Nx-1 else 0
                    jj = j+1 if j < Ny-1 else 0
                    kk = k+1 if k < Nz-1 else 0

                Jp1[0] = Jx[ii,jj,kk]                
                Jp1[1] = Jy[ii,jj,kk]                
                Jp1[2] = Jz[ii,jj,kk]                

                JYx[i,j,k] = (J0[0] + Jp1[0])/2.0
                JYy[i,j,k] = (J0[1] + Jp1[1])/2.0
                JYz[i,j,k] = (J0[2] + Jp1[2])/2.0


    return


def filter_current(W, axis):

    Jp1 = np.zeros(3)
    Jm1 = np.zeros(3)
    J0 = np.zeros(3)
    
    xsweep = False
    ysweep = False
    zsweep = False

    if axis == 0:
        xsweep = True
    if axis == 1:
        ysweep = True
    if axis == 2:
        zsweep = True
     

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                J0[0] = JYx[i,j,k]
                J0[1] = JYy[i,j,k]
                J0[2] = JYz[i,j,k]

                #Boundary  conditions and wrapping
                # TODO wrapping in every dimension done automatically here
                # correct this to check for Nx/y/z_wrap
                ii = i
                jj = j
                kk = k
                if xsweep:
                    ii = xbc(i+1)
                if ysweep:
                    jj = ybc(j+1)

                #if oneD:
                #    ii = i+1 if i < Nx-1 else 0
                #    jj = j
                #    kk = k
                #elif twoD:
                #    ii = i+1 if i < Nx-1 else 0
                #    jj = j+1 if j < Ny-1 else 0
                #    kk = k
                #elif threeD:
                #    ii = i+1 if i < Nx-1 else 0
                #    jj = j+1 if j < Ny-1 else 0
                #    kk = k+1 if k < Nz-1 else 0


                Jp1[0] = JYx[ii,jj,kk]                
                Jp1[1] = JYy[ii,jj,kk]                
                Jp1[2] = JYz[ii,jj,kk]                

                #now get the latter i-1 neighbor
                iii = i
                jjj = j
                kkk = k
                if xsweep:
                    iii = xbc(i-1)
                if ysweep:
                    jjj = ybc(j-1)

                Jm1[0] = JYx[iii,jjj,kkk]                
                Jm1[1] = JYy[iii,jjj,kkk]                
                Jm1[2] = JYz[iii,jjj,kkk]                

                #W = 0.5
                JYx[i,j,k] = (W*Jm1[0] + J0[0] + W*Jp1[0])/(1+2*W)
                JYy[i,j,k] = (W*Jm1[1] + J0[1] + W*Jp1[1])/(1+2*W)
                JYz[i,j,k] = (W*Jm1[2] + J0[2] + W*Jp1[2])/(1+2*W)


    return


#create similar small neighboring cubes as in c++ version
fieldsEx = np.zeros(27)
fieldsEy = np.zeros(27)
fieldsEz = np.zeros(27)

fieldsBx = np.zeros(27)
fieldsBy = np.zeros(27)
fieldsBz = np.zeros(27)

def EB_cube(i,j,k):

    ijk = 0
    for zdir in [-1,0,1]:
        kk = k + zdir 
        if kk < 0:
            kk += Nz
        if kk >= Nz:
            kk -= Nz

        for ydir in [-1, 0, 1]:
            jj = j + ydir 
            if jj < 0:
                jj += Ny
            if jj >= Ny:
                jj -= Ny

            for xdir in [-1, 0, 1]:
                ii = i + xdir 
                if ii < 0:
                    ii += Nx
                if ii >= Nx:
                    ii -= Nx

                if oneD:
                    fieldsEx[ijk] = Ex[ii,j,k]

                    fieldsBx[ijk] = Bx[ii,j,k]
                if twoD:
                    fieldsEx[ijk] = Ex[ii,jj,k]
                    fieldsEy[ijk] = Ey[ii,jj,k]

                    fieldsBx[ijk] = Bx[ii,jj,k]
                    fieldsBy[ijk] = By[ii,jj,k]
                if threeD:
                    fieldsEx[ijk] = Ex[ii,jj,kk]
                    fieldsEy[ijk] = Ey[ii,jj,kk]
                    fieldsEz[ijk] = Ez[ii,jj,kk]

                    fieldsBx[ijk] = Bx[ii,jj,kk]
                    fieldsBy[ijk] = By[ii,jj,kk]
                    fieldsBz[ijk] = Bz[ii,jj,kk]
                ijk += 1

    return

#wrapper functions for global fieldsE and fieldsB

#ijk gives the (Morton) z-ordering of von Neumann neighbors
def ijk(i,j,k):
    return (i+1) + (j+1)*3 + (k+1)*3*3

def exY(i,j,k):
    return fieldsEx[ijk(i,j,k)]
def eyY(i,j,k):
    return fieldsEy[ijk(i,j,k)]
def ezY(i,j,k):
    return fieldsEz[ijk(i,j,k)]

def bxY(i,j,k):
    return fieldsBx[ijk(i,j,k)]
def byY(i,j,k):
    return fieldsBy[ijk(i,j,k)]
def bzY(i,j,k):
    return fieldsBz[ijk(i,j,k)]


# First order staggered grid interpolation
def trilinear_staggered(xd,yd,zd):
    
    # from nodal points
    # f(i+dx) = f(i) + dx * (f(i+1)-f(i))
    # to staggered grid stored at midpoints
    # f at location i+dx  = half of f(i)+f(i-1) + dx*(f(i+1)-f(i-1))
    # where now f(i) means f at location i+1/2.
    # Then we apply the normal linear volume interpolation
    # Note that E and B differ in staggered grid locations

    ex = ( (1 - zd)*((1 - yd)*((0.5 - 0.5*xd)*exY(-1,0,0) + 0.5*exY(0,0,0) + 0.5*xd*exY(1,0,0)) + 
            yd*((0.5 - 0.5*xd)*exY(-1,1,0) + 0.5*exY(0,1,0) + 0.5*xd*exY(1,1,0))) + 
        zd*((1 - yd)*((0.5 - 0.5*xd)*exY(-1,0,1) + 0.5*exY(0,0,1) + 0.5*xd*exY(1,0,1)) + 
                yd*((0.5 - 0.5*xd)*exY(-1,1,1) + 0.5*exY(0,1,1) + 0.5*xd*exY(1,1,1))) )

    ey = ( (1 - zd)*((1 - yd)*(-0.5*(-1. + xd)*(eyY(0,-1,0) + eyY(0,0,0)) + 0.5*xd*(eyY(1,-1,0) + eyY(1,0,0))) + 
            yd*(-0.5*(-1. + xd)*(eyY(0,0,0) + eyY(0,1,0)) + 0.5*xd*(eyY(1,0,0) + eyY(1,1,0)))) + 
        zd*((1 - yd)*(-0.5*(-1. + xd)*(eyY(0,-1,1) + eyY(0,0,1)) + 0.5*xd*(eyY(1,-1,1) + eyY(1,0,1))) + 
                yd*(-0.5*(-1. + xd)*(eyY(0,0,1) + eyY(0,1,1)) + 0.5*xd*(eyY(1,0,1) + eyY(1,1,1)))) )

    ez = ( (1 - zd)*((1 - yd)*(-0.5*(-1. + xd)*(ezY(0,0,-1) + ezY(0,0,0)) + 0.5*xd*(ezY(1,0,-1) + ezY(1,0,0))) + 
            yd*(-0.5*(-1. + xd)*(ezY(0,1,-1) + ezY(0,1,0)) + 0.5*xd*(ezY(1,1,-1) + ezY(1,1,0)))) + 
        zd*((1 - yd)*(-0.5*(-1. + xd)*(ezY(0,0,0) + ezY(0,0,1)) + 0.5*xd*(ezY(1,0,0) + ezY(1,0,1))) + 
                yd*(-0.5*(-1. + xd)*(ezY(0,1,0) + ezY(0,1,1)) + 0.5*xd*(ezY(1,1,0) + ezY(1,1,1)))) )




    bx = ( (1 - zd)*((1 - yd)*(-0.25*(-1. + xd)*(bxY(0,-1,-1) + bxY(0,-1,0) + bxY(0,0,-1) + bxY(0,0,0)) + 
                0.25*xd*(bxY(1,-1,-1) + bxY(1,-1,0) + bxY(1,0,-1) + bxY(1,0,0))) + 
            yd*(-0.25*(-1. + xd)*(bxY(0,0,-1) + bxY(0,0,0) + bxY(0,1,-1) + bxY(0,1,0)) + 
                0.25*xd*(bxY(1,0,-1) + bxY(1,0,0) + bxY(1,1,-1) + bxY(1,1,0)))) + 
        zd*((1 - yd)*(-0.25*(-1. + xd)*(bxY(0,-1,0) + bxY(0,-1,1) + bxY(0,0,0) + bxY(0,0,1)) + 
                    0.25*xd*(bxY(1,-1,0) + bxY(1,-1,1) + bxY(1,0,0) + bxY(1,0,1))) + 
                yd*(-0.25*(-1. + xd)*(bxY(0,0,0) + bxY(0,0,1) + bxY(0,1,0) + bxY(0,1,1)) + 
                    0.25*xd*(bxY(1,0,0) + bxY(1,0,1) + bxY(1,1,0) + bxY(1,1,1)))) )

    by = ( (1 - zd)*((1 - yd)*(-0.25*(-1. + xd)*(byY(-1,0,-1) + byY(-1,0,0) + byY(0,0,-1) + byY(0,0,0)) + 
                0.25*xd*(byY(0,0,-1) + byY(0,0,0) + byY(1,0,-1) + byY(1,0,0))) + 
            yd*(-0.25*(-1. + xd)*(byY(-1,1,-1) + byY(-1,1,0) + byY(0,1,-1) + byY(0,1,0)) + 
                0.25*xd*(byY(0,1,-1) + byY(0,1,0) + byY(1,1,-1) + byY(1,1,0)))) + 
        zd*((1 - yd)*(-0.25*(-1. + xd)*(byY(-1,0,0) + byY(-1,0,1) + byY(0,0,0) + byY(0,0,1)) + 
                    0.25*xd*(byY(0,0,0) + byY(0,0,1) + byY(1,0,0) + byY(1,0,1))) + 
                yd*(-0.25*(-1. + xd)*(byY(-1,1,0) + byY(-1,1,1) + byY(0,1,0) + byY(0,1,1)) + 
                    0.25*xd*(byY(0,1,0) + byY(0,1,1) + byY(1,1,0) + byY(1,1,1)))) )

    bz = ( (1 - zd)*((1 - yd)*(-0.25*(-1. + xd)*(bzY(-1,-1,0) + bzY(-1,0,0) + bzY(0,-1,0) + bzY(0,0,0)) + 
                0.25*xd*(bzY(0,-1,0) + bzY(0,0,0) + bzY(1,-1,0) + bzY(1,0,0))) + 
            yd*(-0.25*(-1. + xd)*(bzY(-1,0,0) + bzY(-1,1,0) + bzY(0,0,0) + bzY(0,1,0)) + 
                0.25*xd*(bzY(0,0,0) + bzY(0,1,0) + bzY(1,0,0) + bzY(1,1,0)))) + 
        zd*((1 - yd)*(-0.25*(-1. + xd)*(bzY(-1,-1,1) + bzY(-1,0,1) + bzY(0,-1,1) + bzY(0,0,1)) + 
                    0.25*xd*(bzY(0,-1,1) + bzY(0,0,1) + bzY(1,-1,1) + bzY(1,0,1))) + 
                yd*(-0.25*(-1. + xd)*(bzY(-1,0,1) + bzY(-1,1,1) + bzY(0,0,1) + bzY(0,1,1)) + 
                    0.25*xd*(bzY(0,0,1) + bzY(0,1,1) + bzY(1,0,1) + bzY(1,1,1)))) )


    return ex, ey, ez, bx, by, bz


#Boris pusher and particle propagator
# note: we update every particle at once in vector format
def update_velocities(grid):

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                cell = grid[i,j,k]
                particles = cell.particles

                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                #q = 1.0
                m = me

                x = particles[:,0]
                y = particles[:,1]
                z = particles[:,2]
                ux = particles[:,3]
                uy = particles[:,4]
                uz = particles[:,5]

                qmw = particles[:,6]

                xd = (x - xmin)/dx
                yd = (y - xmin)/dy
                zd = (z - xmin)/dz

                #interpolate E and B
                if zeroD:
                    Exi = Ex[i,j,k]
                    Eyi = Ey[i,j,k]
                    Ezi = Ez[i,j,k]
                    Bxi = Bx[i,j,k]
                    Byi = By[i,j,k]
                    Bzi = Bz[i,j,k]
                else:
                    EB_cube(i,j,k)
                    Exi, Eyi, Ezi, Bxi, Byi, Bzi = trilinear_staggered(xd,yd,zd)


                #old independent update
                #uxm = ux + q*e*Exi*dt/(2.0*m*c)
                #uym = uy + q*e*Eyi*dt/(2.0*m*c)
                #uzm = uz + q*e*Ezi*dt/(2.0*m*c)

                #update with weights
                # XXX divide by m or not?
                uxm = ux + qmw*Exi*dt/(2.0*m*c)
                uym = uy + qmw*Eyi*dt/(2.0*m*c)
                uzm = uz + qmw*Ezi*dt/(2.0*m*c)

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



#Vay pusher and particle propagator
def Vay_update_velocities(grid):

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                cell = grid[i,j,k]
                particles = cell.particles

                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                #FIXME assuming electron-positron pair plasma
                #q = 1.0
                m = me

                x = particles[:,0]
                y = particles[:,1]
                z = particles[:,2]
                u = particles[:,3]
                v = particles[:,4]
                w = particles[:,5]

                qwe = particles[:,6]

                xd = (x - xmin)/dx
                yd = (y - xmin)/dy
                zd = (z - xmin)/dz

                #interpolate E and B
                if zeroD:
                    Exi = Ex[i,j,k]
                    Eyi = Ey[i,j,k]
                    Ezi = Ez[i,j,k]
                    Bxi = Bx[i,j,k]
                    Byi = By[i,j,k]
                    Bzi = Bz[i,j,k]
                else:
                    EB_cube(i,j,k)
                    Exi, Eyi, Ezi, Bxi, Byi, Bzi = trilinear_staggered(xd,yd,zd)

                #Calculate auxiliary quantities
                cinv = 1.0/c
                #qm = qe/me #electron charge/mass
                q = np.sign(qwe) #extract sign and dispose numerical weight
                qm = q/me #charge to mass ratio assuming pair plasma

                #sign is taken into account only later on
                Bxi *= 0.5*cinv
                Byi *= 0.5*cinv
                Bzi *= 0.5*cinv

                Exi *= 0.5
                Eyi *= 0.5
                Ezi *= 0.5


                #add external fields
                bx0 = Bxi + Bx_ext[i,j,k]*0.5*cinv
                by0 = Byi + By_ext[i,j,k]*0.5*cinv
                bz0 = Bzi + Bz_ext[i,j,k]*0.5*cinv

                ex0 = Exi + Ex_ext[i,j,k]*0.5
                ey0 = Eyi + Ey_ext[i,j,k]*0.5
                ez0 = Ezi + Ez_ext[i,j,k]*0.5



                #First half
                gamma = 1.0/sqrt(1.0 + u*u + v*v + w*w) #reciprocal of the Lorentz force
                #gamma = 1.0

                vx0 = c*u*gamma #3-velocity
                vy0 = c*v*gamma
                vz0 = c*w*gamma

                #u', where cinv is already incorporated
                u1 = c*u + 2.0*ex0*qm + vy0*bz0*qm - vz0*by0*qm
                v1 = c*v + 2.0*ey0*qm + vz0*bx0*qm - vx0*bz0*qm
                w1 = c*w + 2.0*ez0*qm + vx0*by0*qm - vy0*bx0*qm

                #Lorentz factor for u'
                ustar = cinv*(u1*bx0*qm + v1*by0*qm + w1*bz0*qm)
                sigma = cinv*cinv*(c**2 + u1*u1 + v1*v1 + w1*w1) - (bx0*bx0 + by0*by0 + bz0*bz0)
                gamma = 1.0/sqrt(0.5 * (sigma + sqrt(sigma*sigma + 4.0*(bx0*bx0 + by0*by0 + bz0*bz0 + ustar*ustar))))
                #gamma = 1.0
                
                tx = bx0*qm*gamma
                ty = by0*qm*gamma
                tz = bz0*qm*gamma
                f = 1.0/(1.0 + tx*tx + ty*ty + tz*tz)

                u0 = f*(u1 + (u1*tx + v1*ty + w1*tz)*tx + v1*tz - w1*ty)
                v0 = f*(v1 + (u1*tx + v1*ty + w1*tz)*ty + w1*tx - u1*tz)
                w0 = f*(w1 + (u1*tx + v1*ty + w1*tz)*tz + u1*ty - v1*tx)


                #Get normalized 4-velocity
                particles[:,3] = u0*cinv
                particles[:,4] = v0*cinv
                particles[:,4] = w0*cinv


                #Advance position (this is in fact c/gamma)
                # where the extra c takes into account that u=u/c in the position update
                gamma = c*c/sqrt(c*c + u0*u0 + v0*v0 + w0*w0)
                #gamma = c

                if oneD:
                    particles[:,0] += particles[:,3]*gamma*dt
                if twoD:                                                    
                    particles[:,0] += particles[:,3]*gamma*dt
                    particles[:,1] += particles[:,4]*gamma*dt
                if threeD:                                                 
                    particles[:,0] += particles[:,3]*gamma*dt
                    particles[:,1] += particles[:,4]*gamma*dt
                    particles[:,2] += particles[:,5]*gamma*dt




    return



# sort particles between neighboring grid cells
def sort_particles_between_cells(grid):

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):

                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                ip = 0
                while (ip < grid[i,j,k].Npe):
                    particle = grid[i,j,k].particles[ip,:]
                    iddf = (particle[0] - xmin - dx)/dx
                    jddf = (particle[1] - ymin - dy)/dy
                    kddf = (particle[2] - zmin - dz)/dz

                    idd = int(iddf) if iddf < 0 else int(iddf+1)
                    jdd = int(jddf) if jddf < 0 else int(jddf+1)
                    kdd = int(kddf) if kddf < 0 else int(kddf+1)

                    if idd == 0 and jdd == 0 and kdd == 0:
                        #print "inside"
                        ip += 1
                    else:
                        newi = i + idd
                        newj = j + jdd
                        newk = k + kdd

                        #periodic boundary conditions
                        remove_particle = False

                        #if newi < 0:
                        while newi < 0:
                            if Nx_wrap:
                                newi += Nx
                                particle[0] += grid_xmax
                            else:
                                remove_particle = True
                                break
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
                            grid[newi, newj, newk].particles = np.concatenate( (grid[newi,newj,newk].particles, [particle] ), axis=0)
                            grid[newi, newj, newk].Npe += 1

                        #remove from previous cell
                        grid[i,j,k].particles = np.delete(grid[i,j,k].particles, ip, 0)
                        grid[i,j,k].Npe -= 1
                        ip += 1

    return 


def push_half_B():
    ds = np.zeros(3)
    EY = np.zeros(3)
    BY = np.zeros(3)
    EYp1 = np.zeros(3)

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):

                dx, dy, dz = grid_lengths(i,j,k)
                ds[0] = 1.0/dx
                ds[1] = 1.0/dy
                ds[2] = 1.0/dz

                EY[0] = Ex[i,j,k]                
                EY[1] = Ey[i,j,k]                
                EY[2] = Ez[i,j,k]                
                
                #Boundary  conditions and wrapping
                # TODO wrapping in every dimension done automatically here
                # correct this to check for Nx/y/z_wrap
                # EYp1[:] = 0.0

                if oneD:
                    ii = i+1 if i < Nx-1 else 0
                    jj = j
                    kk = k
                elif twoD:
                    ii = i+1 if i < Nx-1 else 0
                    jj = j+1 if j < Ny-1 else 0
                    kk = k
                elif threeD:
                    ii = i+1 if i < Nx-1 else 0
                    jj = j+1 if j < Ny-1 else 0
                    kk = k+1 if k < Nz-1 else 0

                EYp1[0] = Ex[ii,jj,kk]                
                EYp1[1] = Ey[ii,jj,kk]                
                EYp1[2] = Ez[ii,jj,kk]                

                #BY = (c*dt/2.0) * np.cross(ds, EY - EYp1)
                #we update this as b = c B
                BY = (c*dt/2.0) * np.cross(ds, EY - EYp1)

                #add to the field (and transform B=b/c)
                Bx[i,j,k] -= BY[0]/c
                By[i,j,k] -= BY[1]/c
                Bz[i,j,k] -= BY[2]/c


    return


def push_E():

    ds = np.zeros(3)
    JY = np.zeros(3)
    EY = np.zeros(3)
    BY = np.zeros(3)
    BYm1 = np.zeros(3)

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                dx, dy, dz = grid_lengths(i,j,k)
                ds[0] = 1.0/dx
                ds[1] = 1.0/dy
                ds[2] = 1.0/dz

                BY[0] = Bx[i,j,k]                
                BY[1] = By[i,j,k]                
                BY[2] = Bz[i,j,k]                
                
                #Boundary  conditions and wrapping
                # TODO wrapping in every dimension done automatically here
                # correct this to check for Nx/y/z_wrap
                # EYp1[:] = 0.0

                if oneD:
                    ii = i-1 if i > 0 else Nx-1
                    jj = j
                    kk = k
                elif twoD:
                    ii = i-1 if i < 0 else Nx-1
                    jj = j-1 if j < 0 else Ny-1
                    kk = k
                elif threeD:
                    ii = i-1 if i < 0 else Nx-1
                    jj = j-1 if j < 0 else Ny-1
                    kk = k-1 if k < 0 else Nz-1

                BYm1[0] = Bx[ii,jj,kk]                
                BYm1[1] = By[ii,jj,kk]                
                BYm1[2] = Bz[ii,jj,kk]                

                JY[0] = JYx[i,j,k]
                JY[1] = JYy[i,j,k]
                JY[2] = JYz[i,j,k]

                #E_n+1 = E_n + dt*[ curl B - 4pi J ]
                #EY = (c*dt) * np.cross(ds, BY - BYm1) - 4.0*pi*dt*JY
                EY = (c*c*dt) * (np.cross(ds, BY - BYm1))  - dt*JY

                #add to the field
                Ex[i,j,k] += EY[0]
                Ey[i,j,k] += EY[1]
                Ez[i,j,k] += EY[2]


    return

##################################################

def collect_grid(grid):
    particles = np.empty((0,7), dtype=float64)

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                cell = grid[i,j,k]
                particles = np.concatenate((particles, cell.particles), axis=0)
    return particles

def divide_species(particles):
    electrons = particles[np.where(particles[:,6] > 0)]
    positrons = particles[np.where(particles[:,6] < 0)]

    return electrons, positrons



##################################################
##################################################
##################################################
def test_interp(grid):
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                cell = grid[i,j,k]

                xmin,xmax, ymin,ymax, zmin,zmax = grid_limits(i,j,k)
                dx,dy,dz = grid_lengths(i,j,k)

                print "grid limits"
                print xmin, xmax, dx
                print ymin, ymax, dy
                print zmin, zmax, dz

                #interpolate E and B
                EB_cube(i,j,k)

                y = ymin+dy
                z = zmin+dz
                for x in np.linspace(xmin, xmax, 10):
                    print "  x = ", x, "  |", y, z
                    xd = (x - xmin)/dx
                    yd = (y - xmin)/dy
                    zd = (z - xmin)/dz
                    print "  xd = ", xd,"  |", yd, zd

                    Exi, Eyi, Ezi, Bxi, Byi, Bzi = trilinear_staggered(xd,yd,zd)
                    print "   Ex=", Exi, "  |", Eyi, Ezi
                    print "   Bx=", Bxi, "  |", Byi, Bzi

    return
    


