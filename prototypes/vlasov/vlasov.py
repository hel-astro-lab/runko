import numpy as np
from initial import initial
from charge import charge
from poisson import poisson

from position import position_linear
from velocity import velocity_linear
from current import current
from efield import efield



#initialize
#-------------------------------------------------- 
#load configuration
import twostream as prm

#hdiag = diagnostics_init(prm)



#initial(prm)
#ff, gx, gv, ex, fex, ajx, xx, vx, kx, kv, ifdiag = initial(prm)
ff, gx, gv, ex, fex, ajx, xx, vx, kx, kv = initial(prm)

#initial step
rho = charge(ff, prm)
ex,fex = poisson(ex, rho, prm)
ff = position_linear(ff, vx, prm)
ajx = current(ff, vx, prm)


#conservative schemes 
# 0 linear
# 1 2nd order
# 2 4th order
# 6 CPIC4
ex, fex = efield(ex, ajx, prm)

#non-conservative
# 3 - cubic spline
# 5 - CIP3 
#rho = charge(ff, prm)
#ex, fex = poisson(ex, rho, prm)

#3rd possible case ???
#fex = efield_f(fex, ajx, prm)

#-------------------------------------------------- 
# visualization
from pylab import *
from matplotlib import cm
import os, sys




def plot_phasespace(ax, xxi, vxi, ffi, kk):

    #pick kth species
    xx = xxi
    vx = vxi[:,kk]
    ff = ffi[:,:,kk]

    # xx vx[:, kk], ff[:,:,kk] ff[vx:, x:, species]
    ax.set_xlim(xx[0], xx[-1])
    ax.set_ylim(vx[0], vx[-1])
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$v_{x}$')
    
    X, Y = np.meshgrid(xx, vx)
    ax.pcolor(X, Y, ff, 
            cmap='Reds', 
            vmin=ff.min(),
            vmax=ff.max(),
            )

def plot_mean_velocity_pdf(ax, vx, ff, kk):

    ax.set_xlim(-6, 6)
    #ax.set_ylim(

    ax.set_xlabel(r'$v_{x}$')
    #ax.set_ylabel(r'pdf')

    fv = np.mean(ff[:,3:prm.nx+3, kk], 1)
    ax.plot(vx, fv, "k-")


def plot_field(ax, xx, f, quantity):
    ax.set_xlim(xx[0], xx[-1])
    #ax.set_ylim(vx[0], vx[-1])
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(quantity)

    ax.plot(xx, f, "b-")


def visualize(step):

    path = "out"
    if not os.path.exists(path):
        os.makedirs(path)
    
    fname = path+'/vlasov'+str(step)+'.png'
    
    #set up figure
    fig = figure(figsize=(10, 12))
    rc('font', family='serif')
    rc('xtick', labelsize='xx-small')
    rc('ytick', labelsize='xx-small')
    
    gs = GridSpec(3, 3)
    gs.update(hspace = 0.2)
    gs.update(wspace = 0.2)
    
    ax1a = subplot(gs[0, 0:2])
    ax1b = subplot(gs[0, 2])
    
    ax2a = subplot(gs[1, 0:2])
    ax2b = subplot(gs[1, 2])
    
    ax3a = subplot(gs[2, 0])
    ax3b = subplot(gs[2, 1])
    ax3c = subplot(gs[2, 2])

    #phase spaces
    plot_phasespace(ax1a, xx, vx, ff, 0)
    plot_mean_velocity_pdf(ax1b, vx, ff, 0)
    
    plot_phasespace(ax2a, xx, vx, ff, 1)
    plot_mean_velocity_pdf(ax2b, vx, ff, 1)
    
    #plot fields
    plot_field(ax3a, xx, fex, "fex")
    plot_field(ax3b, xx, rho, "rho")
    plot_field(ax3c, xx, ajx, "ajx")
    
    savefig(fname)




#-------------------------------------------------- 
# main loop

jtime = 0
for jtime in range(prm.ntime):
    print "-----------", jtime, "----------"
    visualize(jtime)

    ff  = velocity_linear(ff, fex, prm)
    ff  = position_linear(ff, vx, prm)

    ajx = current(ff, vx, prm)
    rho = charge(ff, prm)

    ex, fex = efield(ex, ajx, prm)
    











