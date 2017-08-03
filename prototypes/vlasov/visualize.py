import numpy as np
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

    fv = np.mean(ff[:,3:-3, kk], 1)
    ax.plot(vx, fv, "k-")


def plot_field(ax, xx, f, quantity):
    ax.set_xlim(xx[0], xx[-1])
    #ax.set_ylim(vx[0], vx[-1])
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(quantity)

    ax.plot(xx, f, "b-")


def visualize(fig, step, xx, vx, ff, ex, rho, ajx):

    path = "out"
    if not os.path.exists(path):
        os.makedirs(path)
    fname = path+'/vlasov'+str(step)+'.png'
    
    gs = GridSpec(3, 3)
    gs.update(hspace = 0.2)
    gs.update(wspace = 0.4)
    
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
    plot_field(ax3a, xx, ex, r'$E_x$')
    plot_field(ax3b, xx, rho, r'$\rho_q$')
    plot_field(ax3c, xx, ajx, r'$J_x$')
    
    savefig(fname)



