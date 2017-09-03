import numpy as np
import os, sys
import copy
import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm

# Copy & append
def cappend(arr, x):
    tmpa = copy.deepcopy( arr )
    tmpx = copy.deepcopy( x )
    tmpa.append( tmpx )
    return tmpa

# copy delete
def cdel(arr, i):
    tmp = copy.deepcopy( arr )
    del tmp[i]
    return tmp



#physical "real" distribution to compare against
def physical_vel(x,y,z):

    mux = 2.0
    muy = 0.0
    muz = 0.0
    sigmax = 5.0
    sigmay = 5.0
    sigmaz = 5.0

    vx = np.exp(-(x-mux)**2 / sigmax**2 )
    vy = np.exp(-(y-muy)**2 / sigmay**2 )
    vz = np.exp(-(z-muz)**2 / sigmaz**2 )

    return vx*vy*vz


# velocity block with data and pointers to neighbors
class vBlock:

    data = 0.0 # value of f(vx,vy,vz) at location (x,y,z)
    loc  = (0.0, 0.0, 0.0)  #block location
    dls  = (0.0, 0.0, 0.0)  #block dimensions

    refLevel = 0
    
    def __init__(self, vx,vy,vz, dx,dy,dz):

        self.loc = (vx, vy, vz)
        self.dls = (dx, dy, dz)





class vMesh:

    mesh = []


    def __init__(self, mins, maxs, dvs):


        mins += dvs/2.0
        maxs -= dvs/2.0

        zi = mins[2]
        while zi <= maxs[2]:

            yi = mins[1]
            while yi <= maxs[1]:

                xi = mins[0]
                while xi <= maxs[0]:
                    #print "({},{},{})".format(xi, yi, zi)

                    vb      = vBlock( xi, yi, zi , dvs[0], dvs[1], dvs[2] )
                    vb.data =  physical_vel(xi, yi, zi)

                    #self.mesh.append( vb )
                    self.mesh = cappend(self.mesh, vb)

                    xi += dvs[0]
                yi += dvs[1]
            zi += dvs[2]





# 3D adaptive momentum mesh
class amr:

    meshx = []
    vx    = []

    meshDx = []

    def __init__(self):

        self.vmin = (-10.0, 0.0, 0.0)
        self.vmax = (+10.0, 0.0, 0.0)
        self.nvMin = 10

        #create initial guiding grid
        dx = (self.vmax[0] - self.vmin[0])/(self.nvMin -1)
        x = self.vmin[0]
        for i in range(self.nvMin):
            self.vx.append( x )

            y = physical_vel(x, 0.0, 0.0)
            self.meshx.append( y )

            x += dx

            self.meshDx.append( 0.0 )


    def update_deriv(self):
        for i in range(self.nvMin-1):
            dx = self.vx[i+1] - self.vx[i]
            dy = self.meshx[i+1] - self.meshx[i]

            self.meshDx[i] = dy/dx



    def visualize_grid(self, ax):

        #ax.cla()
        ax.plot(self.vx, self.meshx, "k-", marker='.')


    def visualize_deriv(self, ax):

        self.update_deriv()

        #ax.cla()
        ax.plot(self.vx, self.meshDx, "k-", marker='.')



def plot_center(ax, vb):
    x = vb.loc[0]
    y = vb.loc[1]
    z = vb.loc[2]
    ax.plot(x, y, marker='.', color='black')


def get_vertex(vb):
    x = vb.loc[0]
    y = vb.loc[1]
    z = vb.loc[2]
    dx = vb.dls[0]
    dy = vb.dls[1]
    dz = vb.dls[2]
   
    vrtx = np.zeros((2))
    vrty = np.zeros((2))
    vrtz = np.zeros((2))

    vrtx[0] = ( x - dx/2.0 )
    vrtx[1] = ( x + dx/2.0 )

    vrty[0] = ( y - dy/2.0 )
    vrty[1] = ( y + dy/2.0 )

    vrtz[0] = ( z - dz/2.0 )
    vrtz[1] = ( z + dz/2.0 )

    return vrtx, vrty, vrtz


def get_xy_bounding_curve(vb):
    vrtx, vrty, vrtz = get_vertex(vb)

    vrtxs = np.zeros(4)
    vrtys = np.zeros(4)

    vrtxs[0] = vrtx[0]
    vrtxs[1] = vrtx[1]
    vrtxs[2] = vrtx[1]
    vrtxs[3] = vrtx[0]

    vrtys[0] = vrty[0]
    vrtys[1] = vrty[0]
    vrtys[2] = vrty[1]
    vrtys[3] = vrty[1]

    return vrtxs, vrtys


def plot_edges(ax, vb):
    vrtxs, vrtys = get_xy_bounding_curve(vb)

    if vb.refLevel == 0:
        col = 'black'

    ax.plot(vrtxs, vrtys, linestyle='solid', color=col)



def visualize_mesh(ax, mesh):

    ax.cla()
    ax.minorticks_on()
    ax.set_xlim(-10.0, 10.0)
    ax.set_ylim(-10.0, 10.0)

    for vb in mesh.mesh:
        plot_center(ax, vb)

        plot_edges(ax, vb)




if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,4))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 2)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )



    ################################################## 
    # set-up grid
    mins = np.array([ -10.0, -10.0, -1.0 ])
    maxs = np.array([  10.0,  10.0,  1.0 ])
    dvs  = np.array([  4.0,    4.0,  2.0 ])
    mesh = vMesh( mins, maxs, dvs )


    visualize_mesh(axs[0], mesh)


    plt.savefig("amr.png")


if __name__ == "old":

    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,4))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(2, 1)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )


    xmin = -12.0
    xmax =  12.0
    ymin = 0.0
    ymax = 1.0

    for ax in axs:
        ax.minorticks_on()
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    
    axs[1].set_ylim(-0.3, 0.3)

    
    velgrid = amr()

    velgrid.visualize_grid(axs[0])
    velgrid.visualize_deriv(axs[1])


    plt.savefig("amr.png")
