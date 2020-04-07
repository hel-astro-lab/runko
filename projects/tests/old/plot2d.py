import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os
import matplotlib.ticker as ticker
from scipy.stats import mstats
from scipy.optimize import curve_fit

from visualize import imshow

from configSetup import Configuration
from combine_files import get_file_list
from combine_files import combine_tiles



# trick to make nice colorbars
# see http://joseph-long.com/writing/colorbars/
def colorbar(mappable, 
        loc="right", 
        orientation="vertical", 
        size="5%", 
        pad=0.05, 
        ticklocation='auto'):
        #loc="top", 
        #orientation="horizontal", 
        #size="8%", 
        #pad=0.03, 
        #ticklocation='top'):

    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(loc, size=size, pad=pad)
    return fig.colorbar(mappable, cax=cax, orientation=orientation, ticklocation=ticklocation)



def plot2d(ax, val, title="", vmin=None, vmax=None, cmap="RdBu"):
    
    nx, ny = np.shape(val)
    xmin = 0.0
    ymin = 0.0
    xmax = nx
    ymax = ny

    if vmin == None:
        vmin, vmax = np.min(val), np.max(val)
        vminmax = np.maximum( np.abs(vmin), np.abs(vmax) )
        vmin = -vminmax
        vmax =  vminmax
    elif not(vmin == None) and (vmax == None):
        vmax = np.max(val)

    print("2D: vmin: {} vmax: {} {}".format(vmin, vmax, title))

    im = imshow(ax, val, xmin, xmax, ymin, ymax,
           cmap = cmap,
           vmin = vmin,
           vmax = vmax,
           )

    cb = colorbar(im)

    ax.set_title(title)


    return im, cb



plt.fig = plt.figure(1, figsize=(8,10))

plt.rc('font', family='serif', size=7)
plt.rc('xtick')
plt.rc('ytick')

gs = plt.GridSpec(4, 3)
gs.update(hspace = 0.4)
gs.update(wspace = 0.7)

axs = []
for ai in range(9):
    axs.append( plt.subplot(gs[ai]) )

axs.append( plt.subplot(gs[3,0:3]) )


#--------------------------------------------------
# read tristan output files
#import glob
#files=sorted(glob.glob('output/flds.tot*'))
#prtfiles=sorted(glob.glob('output/prtl.tot*'))
#
## go through the files
#start=0
#end=len(files)
##end=20
#interval=1
##global d
##d=[]

#--------------------------------------------------
# read plasmabox output files

conf = Configuration('config-weibel.ini') 
fdir = "weibel/out/"
fname_F = "field"
fname_A = "analysis"

files_F = get_file_list(fdir, fname_F)
files_A = get_file_list(fdir, fname_A)

ex = []
ey = []
ez = []
bx = []
by = []
bz = []
jx = []
jy = []
jz = []
dens = []

start = 0
end = len(files_F)
interval = 1


for step in range(start,end,interval):
    print("reading {} and {}".format(files_F[step], files_A[step]))

    exi = combine_tiles(files_F[step], "ex", conf)[:,:,0]
    eyi = combine_tiles(files_F[step], "ey", conf)[:,:,0]
    ezi = combine_tiles(files_F[step], "ez", conf)[:,:,0]
                                                         
    jxi = combine_tiles(files_F[step], "jx", conf)[:,:,0]
    jyi = combine_tiles(files_F[step], "jy", conf)[:,:,0]
    jzi = combine_tiles(files_F[step], "jz", conf)[:,:,0]
                                                         
    bxi = combine_tiles(files_F[step], "bx", conf)[:,:,0]
    byi = combine_tiles(files_F[step], "by", conf)[:,:,0]
    bzi = combine_tiles(files_F[step], "bz", conf)[:,:,0]

    densi = combine_tiles(files_F[step], "rho", conf)[:,:,0]
    
    print("shape ex")
    print( np.shape(exi) )


    xmin = 0.0
    ymin = 0.0
    xmax = 1.0
    ymax = 1.0

    nx, ny = np.shape(densi)
    print("nx={} ny={}".format(nx, ny))


    im0, cb0 = plot2d(axs[0], exi, title=r"$E_x$")
    im1, cb1 = plot2d(axs[1], eyi, title=r"$E_y$")
    im2, cb2 = plot2d(axs[2], ezi, title=r"$E_z$")

    im3, cb3 = plot2d(axs[3], bxi, title=r"$B_x$")
    im4, cb4 = plot2d(axs[4], byi, title=r"$B_y$")
    im5, cb5 = plot2d(axs[5], bzi, title=r"$B_z$")

    im6, cb6 = plot2d(axs[6], jxi, title=r"$J_x$")
    im7, cb7 = plot2d(axs[7], jyi, title=r"$J_y$")
    im8, cb8 = plot2d(axs[8], jzi, title=r"$J_z$")

    im9, cb9 = plot2d(axs[9],densi, title=r"$\rho$", vmin=0.0, cmap="viridis")


    #ymin = ny/2 - 50
    #ymax = ny/2 + 50
    #for ax in axs:
    #    ax.set_ylim(ymin, ymax)


    slap = str(step).rjust(4, '0')
    fname = fdir + 'visz_{}.png'.format(slap)
    plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.97)
    plt.savefig(fname)

    for cb in [cb0, cb1, cb2, cb3, cb4, cb5, cb6, cb7, cb8, cb9]:
        cb.remove()


    ex.append(exi)
    ey.append(eyi)
    ez.append(ezi)
    jx.append(jxi)
    jy.append(jyi)
    jz.append(jzi)
    bx.append(bxi)
    by.append(byi)
    bz.append(bzi)
    dens.append(densi)

    #pxe   = np.squeeze(f1['ue']).flatten()
    #pye   = np.squeeze(f1['ve']).flatten()
    #pze   = np.squeeze(f1['we']).flatten()
    #xe    = np.squeeze(f1['xe']).flatten()
    #xarr = np.arange(len(ex))/40.0


    #res1a = ax1.plot(xe, pxe, "k.", alpha=0.8)
    #res2a = ax2.plot(xarr, dens,"k-", alpha=0.8)
    #res3a = ax3.plot(xarr, jx,  "k-", alpha=0.8)
    #res4a = ax4.plot(xarr, ex,  "k-", alpha=0.8)

    ##res2 = ax3.plot([0,:,0], rpic.By[:,0,0], "b-")


    #fname = path+'/oneD_'+str(step)+'.png'
    #savefig(fname)

    ##clean figure
    #res1a.pop(0).remove()
    #res2a.pop(0).remove()
    #res3a.pop(0).remove()
    #res4a.pop(0).remove()



# transform into numpy array
ex = np.array(ex) 
ey = np.array(ey) 
ez = np.array(ez) 
bx = np.array(bx) 
by = np.array(by) 
bz = np.array(bz) 
jx = np.array(jx) 
jy = np.array(jx) 
jz = np.array(jx) 
dens = np.array(dens) 


