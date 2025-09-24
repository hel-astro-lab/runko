import numpy as np
import h5py as h5
import sys, os

import matplotlib.pyplot as plt
 
from configSetup import Configuration

sys.path.insert(0,'../')
from setup_tests import Configuration_Test
from setup_tests import lap2time

import pytools


if __name__ == "__main__":
    fig = plt.figure(1, figsize=(3.487, 2.65))
    
    plt.rc('font',  family='sans-serif')
    #plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.0)
    gs.update(wspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    #axs.append( plt.subplot(gs[1,0]) )
    
    for ax in axs:
        ax.minorticks_on()

    if len(axs) > 1:
        axs[0].tick_params(which='x', direction="in")
        axs[0].axes.get_xaxis().set_ticks([])

    #axs[0].set_yscale('log')

    axs[0].set_ylabel(r'$\omega/\omega_{pp}$')
    axs[0].set_xlabel(r'$k \lambda_p$')



    #--------------------------------------------------
    # read and plot data

    confs = [
            "tmp.ini",
            ]

    for conf_filename in confs:
        conf = Configuration_Test(conf_filename)

        # restructure path to point to this dir
        fdir = conf.outdir 
        fdir = fdir.split('/')[1] + '/'
        print(fdir)

        #--------------------------------------------------
        # normalization and units
        qe = np.abs(conf.qe)
        me = np.abs(conf.me)*qe
        c  = conf.cfl 
        print("me: {} qe: {} c: {}".format(me, qe, c))

        Lx, Ly, Lz = conf.NxMesh*conf.Nx, conf.NyMesh*conf.Ny, conf.NzMesh*conf.Nz
        norm_prtcl = Lx*Ly*Lz*conf.ppc*2.0 #pair-plasma so multiply by 2
        norm_prtcl *= conf.gamma*me*c**2 #kinetic energy in flow

        norm_flds = 1.0
        norm_flds /= 2.0*(1.0/conf.stride)*(1.0/conf.stride) #divide by volume of simulation domain

        #3D switch
        if Lz > 1:
            norm_flds /= (1.0/conf.stride)

        # fields into units of particle enthalphy density
        norm_flds /= norm_prtcl

        #--------------------------------------------------
        # read data

        #flds_E = []
        #flds_B = []
        #flds_Bpar = []
        #flds_time = []

        et = np.zeros((conf.Nx*conf.NxMesh, conf.Nt))
        break_lap = 2000

        for i,lap in enumerate(range(0, conf.Nt+1, conf.interval)):
            info = {}
            info['lap'] = lap
            info['fields_file'  ] = fdir + 'flds_'+str(lap)+'.h5'

            #print(info['fields_file'])
            if not(os.path.isfile(info['fields_file'])):
                continue

            if lap > break_lap:
                break

            #print("lap {}".format(lap))
            time = lap2time(lap, conf) #, do_print=True)

            # read file
            f5_fld = h5.File(info['fields_file'],'r')
            #bx  = pytools.read_h5_array(f5_fld, 'bx')
            #by  = pytools.read_h5_array(f5_fld, 'by')
            #bz  = pytools.read_h5_array(f5_fld, 'bz')

            #jx  = pytools.read_h5_array(f5_fld, "jx")
            #jy  = pytools.read_h5_array(f5_fld, "jy")
            #jz  = pytools.read_h5_array(f5_fld, "jz")

            ex  = pytools.read_h5_array(f5_fld, "ex")
            #ey  = pytools.read_h5_array(f5_fld, "ey")
            #ez  = pytools.read_h5_array(f5_fld, "ez")

            #print(np.shape(ex))

            et[:,i] = ex[:,0,0] #take first slice

            #end of this time step; close f5
            f5_fld.close()


        x = np.arange(0, conf.Nx*conf.NxMesh)/conf.c_omp
        y = lap2time(np.arange(0, conf.Nt+1, conf.interval), conf)
        y *= np.sqrt(2)

        #y /= 4.0

        y = y[:break_lap]
        et = et[:,:break_lap]

        #axs[0].plot(x, et[:,0], 'r',)
        #axs[0].plot(x, et[:,50], 'b',)
        #axs[0].plot(x, et[:,90], 'g',)

        dt = y[1]-y[0]
        dx = x[1] - x[0]
        print("dx: ", dx)
        print("dt: ", dt)

        #et = et[:, 25:225]
    
        #--------------------------------------------------
        # filter data with windowing
        A = et.T

        if False:
            XX, YY = np.meshgrid(x,y)
            print("shape XX", np.shape(XX))
            print("shape YY", np.shape(YY))
            Lx = conf.Nx*conf.NxMesh
            Nt = conf.Nt/conf.interval
            print("lx", Lx)
            print("nt", Nt)
    
            kx = 2.0*np.pi/Lx
            wy = 2.0*np.pi/Nt

            #A[:,:] = np.sin(kx*XX*100)
            A[:,:] = 0.0
            #for m in range(1,20):
            for m in [150]:
                #A[:,:] += np.sin(wy*YY*m + 2.0*np.pi*np.random.rand())

                A[:,:] += np.sin(kx*XX*m + 2.0*np.pi*np.random.rand())

                #A[:,:] += np.sin(wy*YY*m + 2.0*np.pi*np.random.rand())

                #A[:,:] += m*np.sin(wy*YY*m + 2.0*np.pi*np.random.rand())

                #A[:,:] += np.sin(wy*YY*m)


        (lines, cols) = np.shape(A)
        #window = np.hamming(lines).reshape(lines, 1)
        #window = np.bartlett(lines).reshape(lines, 1)
        #window = np.blackman(lines).reshape(lines, 1)
        window = np.kaiser(lines, 14.0).reshape(lines, 1)
        A *= window
        
        #--------------------------------------------------
        # Fourier transform in space and time
        A_fft = np.fft.rfft2(A)

        (ny, nx) = np.shape(A_fft)
        print("shape after transform:", np.shape(A_fft))

        # spatial wave vector x component k_x (only half is considered due to Nqyust frequency cut)
        dk = 2*np.pi/(nx * dx)

        #dk = 1.0
        k = np.arange(nx)*dk
        k1 = 1
        k2 = nx

        # temporal angular frequency \omega
        dw = 2*np.pi/(ny*dt)

        #dw = 1.0
        w = np.arange(ny)*dw*2
        w1 = 1
        w2 = ny//2


        #compute minimum and maximum values
        Lx = dx*nx
        T  = ny*dt
        
        kmin = 2*np.pi/Lx
        kmax = np.pi/dx
        
        wmin = 2*np.pi/T
        wmax = np.pi/dt
        
        print("x min/max:", dx, Lx)
        print("t min/max:", dt, T)
        print("k min/max:", kmin, kmax)
        print("w min/max:", wmin, wmax)

        # Change to spectra by considering |F] 
        F = np.log10( np.abs(A_fft)**2 )
        #F = np.abs(A_fft)
        #F = A.real**2
        #F = np.log10(A.real)
        #F = A.real

        print("min/max:", np.min(F), np.max(F))
        
        #K, W = np.meshgrid(k, w)
        #print("nx=",nx)
        #print("ny=",ny)
        #print(np.shape(K))
        #print(np.shape(W))

        vmin = -1
        vmax = 3.0

        #faster plotting with imshow
        im = axs[0].imshow(
                F[w1:w2, k1:k2],
                extent=[k[k1], k[k2-1], w[w1], w[w2-1]],
                origin='lower',
                aspect='auto',
                interpolation='nearest',
                cmap='plasma',
                vmin=vmin,
                vmax=vmax,
                #vmin=-0,
                #vmax=2.0
                )

        #axs[0].set_ylim((w[w1],w[w2-1]))
        #axs[0].set_xlim((k[k1],k[k2-1]))
        
        axs[0].set_ylim((0.0, 10.0))
        axs[0].set_xlim((0.0, 10.0))


        

        ##################################################
        #analytical relations
        
        wp = 1.0 #plasma frequency XXX why x2? 
        #vth = np.sqrt(1.0e-5) #thermal velocity XXX why /2?
        #vth = 0.25 #thermal velocity XXX why /2?

        vth = np.sqrt(2.0*conf.delgam) #3D rms velocity


        #numerical propagation speed
        def w_num(k):
            v= dx/dt
            return k*v
        axs[0].plot(k, w_num(k), color='black', linestyle='dotted', linewidth=0.8, alpha=0.5)
        
        
        #electron plasma frequency
        def w_ep(k, wp):
            return np.ones(len(k))*wp
        
        axs[0].plot(k, w_ep(k, wp), color='black', linestyle='dotted', linewidth=0.8, alpha=0.5)
        
        #warm plasma (Langmuir waves)
        def w_warmPlasma(k, wp, vth):
            wps = np.ones(len(k))*wp #*2.0
            vths = 4.0*k*vth #/2.0
            return np.sqrt( wps**2 + 3.0*vths**2 )
        
        axs[0].plot(k, w_warmPlasma(k, wp, vth),
                color='black', 
                linestyle='dashed', 
                linewidth=0.8,
                alpha=0.5)
        



    #--------------------------------------------------    
    # close and save fig
    axleft    = 0.18
    axbottom  = 0.16
    axright   = 0.96
    axtop     = 0.80

    pos1 = axs[0].get_position()
    axwidth  = axright - axleft
    axheight = (axtop - axbottom)*0.02
    axpad = 0.01

    cax = fig.add_axes([axleft, axtop + axpad, axwidth, axheight])
    #cax = fig.add_axes([0.12, 0.86, 0.86, 0.03]) 

    plt.colorbar(im, cax=cax, orientation='horizontal', label=r'log$_{10} | \mathcal{F}( E_x ) |$')
    cax.xaxis.set_label_position('top')
    cax.xaxis.set_ticks_position('top')

    #plt.colorbar(im, cax=cax, orientation='horizontal', label=r'log$_{10} | \mathcal{F}( E_x ) |$')

    #cb1 = colorbar.ColorbarBase(
    #        cax,
    #        cmap='plasma',
    #        norm=norm,
    #        orientation='horizontal',
    #        ticklocation='top')
    #cb1.set_label(r'$t$ $(l_0/c)$')

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    fname = 'wk.pdf'
    plt.savefig(fname)




