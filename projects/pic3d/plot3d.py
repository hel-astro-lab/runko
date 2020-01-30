import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os
import matplotlib.ticker as ticker
from matplotlib import cm

#3D
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable

# pytools bindings
import pytools
from problem import Configuration_Problem as Configuration


def read_h5_array(f5, var_name):

    nx = f5['Nx'].value
    ny = f5['Ny'].value
    nz = f5['Nz'].value

    # column-ordered data with image convention (x horizontal, y vertical, ..)
    # i.e., so-called fortran ordering
    val = f5[var_name][:]

    print("reshaping 1D array of {} into multiD with {} {} {}".format(len(val), nx,ny,nz))


    # reshape to python format; from Fortran image to C matrix
    val = np.reshape(val, (nz, ny, nx))
    val = val.ravel(order='F').reshape((nx,ny,nz))

    #val = val.ravel(order='C').reshape((nx,ny,nz))
    #val = val.reshape(-1) #F to C

    return val



if __name__ == "__main__":

    do_print = True
    do_dark = False

    if do_dark:
        fig = plt.figure(1, figsize=(8,4.0), dpi=300)
        
        plt.rc('font', family='serif', size=7)
        plt.rc('xtick')
        plt.rc('ytick')

        plt.style.use('dark_background')
    else:
        fig = plt.figure(1, figsize=(4,2.5), dpi=200)
        plt.rc('font',  family='sans-serif')
        #plt.rc('text',  usetex=True)
        plt.rc('xtick', labelsize=8)
        plt.rc('ytick', labelsize=8)
        plt.rc('axes',  labelsize=8)
    

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.0)
    gs.update(wspace = 0.0)
    
    axs = []
    #axs.append( plt.subplot(gs[0,0]) )
    axs.append( fig.add_subplot(111, projection='3d') )


    #--------------------------------------------------
    # prepare data

    args = pytools.parse_args()
    conf = Configuration(args.conf_filename, do_print=do_print)

    print(conf.outdir)

    fname_fld = "flds"
    fname_prtcls = "test-prtcls"
    
    #args.lap
    lap = 80

    fields_file   = conf.outdir + '/'+fname_fld+'_'+str(lap)+'.h5'

    f5 = h5.File(fields_file,'r')
    data = read_h5_array(f5, 'jz')

    #cut reflector out
    data = data[6:,:,:]

    #limit box length
    data = data[0:128,:,:]

    print(np.shape(data))

    nx, ny, nz = np.shape(data)

    #--------------------------------------------------
    # create box

    box = pytools.pybox3d.Box(axs[0])

    aspect_ratio = nx/ny
    print("box aspect ratio: {}".format(aspect_ratio))
    box.dx = aspect_ratio
    box.dy = 1.0
    box.dz = 1.5

    #box.set_data(np.log10(data))
    box.set_data(data)
    box.vmin = -0.05
    box.vmax = +0.05
    cmap = cm.get_cmap('RdBu')


    #surface rendering
    box.draw_left( cmap=cmap)
    box.draw_front(cmap=cmap)
    box.draw_top(  cmap=cmap)
    box.draw_outline()



    #back exploded panels
    if True:
        off_bot    = 1.7
        off_back   = 1.7
        off_left   = 0.7
        off_right  = 0.7

        box.draw_exploded_panels_outline("bottom", off=off_bot)
        box.draw_exploded_panels_outline("left",   off=off_left)
        #box.draw_exploded_panels_outline("right",  off=off_right)
        box.draw_exploded_panels_outline("back",   off=off_back)
        
        box.draw_exploded_bottom(off=off_bot,   cmap=cmap)
        box.draw_exploded_back( off=off_back, cmap=cmap)

        box.draw_exploded_left(  off=off_left,  cmap=cmap)
        #box.draw_exploded_right( off=off_right, cmap=cmap)

    if False:
        #front exploded panels
        off = -1.95 #this puts them in front of the box
        cmap = plt.cm.RdBu
        box.draw_exploded_panels_outline("bottom", off=off)
        box.draw_exploded_panels_outline("left",   off=off)
        box.draw_exploded_panels_outline("right",  off=off)
        
        box.draw_exploded_bottom(off=off, cmap=cmap)
        box.draw_exploded_left(  off=off, cmap=cmap)
        box.draw_exploded_right( off=off, cmap=cmap)


    axs[0].set_axis_off()
    axs[0].view_init(45.0, -110.0)

    #axs[0].view_init(45.0, 100.0)


    if False:
        #colorbars
        m = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
        m.set_array([0.0, 1.0])
        m.set_clim( vmin=0.0, vmax=1.0 )
        cbaxes = fig.add_axes([0.2, 0.91, 0.6, 0.02]) #[left, bottom, width, height],
        cb = plt.colorbar(m, cax = cbaxes, orientation="horizontal", ticklocation="top")  
        fig.text(0.15, 0.91,  r'$n_{\pm}$')

        m = plt.cm.ScalarMappable(cmap=plt.cm.inferno)
        m.set_array([0.0, 1.0])
        m.set_clim( vmin=0.0, vmax=1.0 )
        cbaxes = fig.add_axes([0.2, 0.09, 0.6, 0.02]) #[left, bottom, width, height],
        cb = plt.colorbar(m, cax = cbaxes, orientation="horizontal", ticklocation="top")  
        fig.text(0.15, 0.10,  r'$n_{\nu}$')

        m = plt.cm.ScalarMappable(cmap=plt.cm.RdBu)
        m.set_array([-1.0, 1.0])
        m.set_clim( vmin=-1.0, vmax=1.0 )
        cbaxes = fig.add_axes([0.2, 0.06, 0.6, 0.02]) #[left, bottom, width, height],
        cb = plt.colorbar(m, cax = cbaxes, orientation="horizontal", ticklocation="bottom")  
        #cb.set_label(r'$J$', rotation=0)
        fig.text(0.15, 0.05,  r'$J$')


    pytools.pybox3d.axisEqual3D(axs[0])

    #axs[0].set_title('Step {}'.format(lap+1))

    slap = str(lap).rjust(4, '0')
    fname = conf.outdir + '/3d_'+slap
    plt.subplots_adjust(left=-0.45, bottom=-0.45, right=1.35, top=1.45)

    axs[0].set_xlabel('x')
    axs[0].set_ylabel('y')
    axs[0].set_zlabel('z')

    #plt.savefig(fname+'.pdf')
    plt.savefig(fname+'.png')
    


