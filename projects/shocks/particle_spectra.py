import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os
from scipy.special import kv

from configSetup import Configuration
from combine_files import get_file_list
from combine_files import combine_tiles

from parser import parse_input
import argparse


default_values = {
    'xmin': None,
    'xmax': None,
    'clip': None,
    'isp' : None,
    'nbins': 100
}


default_turbulence_values = {
        'gammam1': {'title': r"$\rho$",
                'xmin': 0.1,
                'xmax': 200.0,
                },
        'gamma': {'title': r"$\gamma$",
                'xmin': 1.0,
                'xmax': 1000.0,
                },
}




def create_spectra(
        var,
        info, 
        conf,
        isp  = None,
        xmin = None,
        xmax = None,
        clip = None,
        read_file = False,
        no_plots = False,
        ):

    print("fdir: ", info['fdir'])
    print("lap: ", info['lap'])

    #--------------------------------------------------
    # unpack incoming arguments that modify defaults
    args = {}

    #general defaults
    for key in default_values:
        args[key] = default_values[key]

    #overwrite with turbulence defaults
    try:
        for key in default_turbulence_values[var]:
            args[key] = default_turbulence_values[var][key]
    except:
        pass

    #finally, overwrite with user given values
    for key in args:
        try:
            user_val = eval(key)
            if not(user_val == None):
                args[key] = user_val
                print("overwriting {} key with {}".format(key, user_val))
        except:
            pass

    if args['xmin'] == None:
        args['xmin'] = 1.0
    if args['xmax'] == None:
        args['xmax'] = 1.0e4

    # finally, re-check that user did not give xmin/xmax
    args['xmin'] = xmin if not(xmin == None) else args['xmin']
    args['xmax'] = xmax if not(xmax == None) else args['xmax']

    #nor the project default
    args['xmin'] = default_turbulence_values[var]['xmin'] if not(xmin == None) else args['xmin']
    args['xmax'] = default_turbulence_values[var]['xmax'] if not(xmax == None) else args['xmax']

    #--------------------------------------------------
    # print all options
    print("--- {} ---".format(var))
    for key in args:
        if not(key == None):
            print(" setting {}: {}".format(key, args[key]))


    #--------------------------------------------------
    # actual processing of files

    print('--------------------------------------------------')
    print("reading file set {}".format(info['particle_file']))


    bins = np.logspace(
            np.log10(args['xmin']),
            np.log10(args['xmax']),
            args['nbins'])

    hist = np.zeros(args['nbins']-1)

    histx = np.zeros(args['nbins']-1)
    histy = np.zeros(args['nbins']-1)
    histz = np.zeros(args['nbins']-1)

    if read_file:
        slap = str(info['lap']).rjust(4, '0')
        f5in = h5.File(info['fdir'] + 'prtcl_spec_{}.h5'.format(slap), 'r')
        bins = f5in['bins'][:]
        hist = f5in['hist'][:]
        f5in.close()
    else:
        #build spectra first
        
        # all simulation particle spectra
        if False:
            rank = 0
            while True:
                fname = info['particle_file']+'-'+str(rank)+'_'+str(info['lap'])+'.h5'
                print("{} ".format(rank), end='', flush=True)
            
                if os.path.exists(fname):
                    f5 = h5.File(fname,'r')
                    rank += 1
                else:
                    print("could not find {}; finishing...".format(fname))
                    break

                for f5_tile in f5:
                    for dset in f5[f5_tile]:

                        #skip if not correct species
                        if not(args['isp']==None):
                            if not(f[f5_tile][dset]['ispcs'].value == args['isp']):
                                continue

                        i = f5[f5_tile][dset]['i'].value
                        j = f5[f5_tile][dset]['j'].value
                        k = f5[f5_tile][dset]['k'].value

                        vx = f5[f5_tile][dset]['vx'][:]
                        vy = f5[f5_tile][dset]['vx'][:]
                        vz = f5[f5_tile][dset]['vz'][:]

                    #print("{} particles in tile ({},{},{})".format(len(vx), i,j,k))
                    #gam = np.sqrt(1.0 + vx*vx + vy*vy + vz*vz) - 1.0
                    gam = np.sqrt(1.0 + vx*vx + vy*vy + vz*vz)

                    tmphist, edges = np.histogram(gam, bins=bins)
                    hist += tmphist

                    tmphistx, edges = np.histogram(ux, bins=bins)
                    tmphisty, edges = np.histogram(uy, bins=bins)
                    tmphistz, edges = np.histogram(uz, bins=bins)

                    histx += tmphistx
                    histy += tmphisty
                    histz += tmphistz

        # test particle spectra
        if True:
            fname = info['particle_file'] + "_" + str(info['lap']) + ".h5"
            f5F = h5.File(fname,'r')
            nx = f5F['Nx'].value
            ny = f5F['Ny'].value
            nz = f5F['Nz'].value
            print("reshaping 1D array into multiD with ({} {} {}) {}".format(nx,ny,nz, info['lap']))

            ux   = np.reshape( f5F['vx'][:]   ,(nx,-1)) 
            uy   = np.reshape( f5F['vy'][:]   ,(nx,-1)) 
            uz   = np.reshape( f5F['vz'][:]   ,(nx,-1)) 
            #wgt  = np.reshape( f5F['wgt'][:]  ,(nx,-1)) 

            #gam = np.sqrt(1.0 + ux*ux + uy*uy + uz*uz) - 1.0
            gam = np.sqrt(1.0 + ux*ux + uy*uy + uz*uz) 
            
            tmphist, edges = np.histogram(gam, bins=bins)
            hist += tmphist
            
            tmphistx, edges = np.histogram(ux, bins=bins)
            tmphisty, edges = np.histogram(uy, bins=bins)
            tmphistz, edges = np.histogram(uz, bins=bins)

            histx += tmphistx
            histy += tmphisty
            histz += tmphistz

        #--------------------------------------------------
        # save to file; flushes after every rank
        slap = str(info['lap']).rjust(4, '0')
        f5out = h5.File(info['fdir'] + 'prtcl_spec_{}.h5'.format(slap), 'w')
        dset1  = f5out.create_dataset("bins",  data=bins)
        dset2  = f5out.create_dataset("hist",  data=hist)

        dset2x = f5out.create_dataset("histx",  data=histx)
        dset2y = f5out.create_dataset("histy",  data=histy)
        dset2z = f5out.create_dataset("histz",  data=histz)

        f5out.close()


    if no_plots:
        return bins[:-1], hist


    #--------------------------------------------------
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')

    axs[0].set_ylim((1.0e-1, 3.0e3))
    axs[0].set_xlim((1.0, 2.0e2))
    #axs[0].set_ylim((1.0e3, 1.0e6))

    axs[0].plot(bins[:-1], hist, color='k')

    #axs[0].plot(bins[:-1], histx, color='b', alpha=0.6)
    #axs[0].plot(bins[:-1], histy, color='g', alpha=0.6)
    #axs[0].plot(bins[:-1], histz, color='r', alpha=0.6)


    #xt, yt, kt = fit_temperature(bins, hist)
    #axs[0].plot(xt, yt, "k", linestyle='dotted')


    #xp, yp = fit_slope(bins, hist)
    #yp *= 2.0
    #axs[0].plot(xp, yp, "k--")


    #--------------------------------------------------
    # styling
    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$dN/d\ln(\gamma-1)$")

    #--------------------------------------------------
    axleft    = 0.18
    axbottom  = 0.18
    axright   = 0.96
    axtop     = 0.96
    plt.fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    fname = info['fdir'] + 'prtcl_spec_{}.pdf'.format(slap)
    plt.savefig(fname)


    return bins[:-1], hist



# fit powerlaw slope
def fit_slope(xx, yy):

    x_med = xx[ np.argmax(yy) ]
    print("median gamma {}".format(x_med))
    x1 = 15.0
    x2 = 30.0

    #x1 = 3.0
    #x2 = 7.0

    #x1 = 8.0
    #x2 = 20.0

    i1 = 0
    for i,x in enumerate(xx):
        if x > x1:
            i1 = i
            break

    i2 = len(xx)
    for i in range(len(xx)-1, 0, -1):
        if xx[i] < x2:
            i2 = i
            break

    print("limits are {} {} / {} {}".format(i1, i2, xx[i1], xx[i2]))

    xx1 = xx[i1:i2]
    yy1 = yy[i1:i2]

    xf = np.log(xx1[np.nonzero(yy1)])
    yf = np.log(yy1[np.nonzero(yy1)])
    
    fit = np.polyfit(xf, yf, 1)
    print("slope of {}".format(fit[0]))

    poly = np.poly1d(fit)

    xx2 = np.exp(xf)
    yy2 = np.exp(poly(xf))

    return xx2, yy2


def maxwell_juttner(gam, kt):
    #gam = gamm1 + 1.0
    #b2 = kv(2, 1.0/kt)
    beta = np.sqrt(1.0 - 1.0/gam**2)
    #p1 = gam*gam*beta/kt/b2
    #p1 = beta*gam**2.0
    f = gam*beta*np.exp(-gam/kt)
    f /= np.max(f)
    
    return f

    #dg = gam
    #return gam*gam*dg*beta*np.exp(-gam/kt)



# fit relativistic Maxwell-Juttner 
def fit_temperature(xx, yy):
    x_med = xx[ np.argmax(yy) ]
    i1 = 0
    x2 = 5.0*x_med

    i2 = len(xx)
    for i in range(len(xx)-1, 0, -1):
        if xx[i] < x2:
            i2 = i
            break
    print("limits are {} {} / {} {}".format(i1, i2, xx[i1], xx[i2]))

    xf = xx[i1:i2]
    #yf = maxwell_juttner(xf, 3.3)
    yf = maxwell_juttner(xf, 0.4)

    yf *= 0.9*np.max(yy)/np.max(yf)


    return xf, yf, 0.1



#--------------------------------------------------

if __name__ == "__main__":


    #plt.fig = plt.figure(1, figsize=(8,7), dpi=300)
    plt.fig = plt.figure(1, figsize=(3.487, 2.2))
    
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
    
    for ax in axs:
        ax.minorticks_on()

    #--------------------------------------------------
    # command line driven version

    conf, fdir, args = parse_input()

    fdir += '/'
    print("plotting {}".format(fdir))
    
    #fname_P = "particles"
    fname_P = "test-prtcls"


    # if lap is defined, only process that one individual round
    if not(args.lap == None):
        info = {}
        info['lap'] = args.lap
        info['fdir'] = fdir
        #info['particle_file'] = fdir+'restart/'+fname_P
        info['particle_file'] = fdir+fname_P

        #create_spectra('gamma', info, conf)
        create_spectra('gamma', info, conf, read_file=False)

    # else process every file that there is
    else:
        files_P = []
        for lap in range(conf.interval, conf.Nt+1, conf.interval):
            files_P.append( fdir+fname_P)
    
            #info = build_info(fdir, lap)
            info = {}
            info['lap'] = lap
            info['fdir'] = fdir
            info['particle_file'] = fdir+fname_P

            create_spectra('gamma', info, conf)

    
