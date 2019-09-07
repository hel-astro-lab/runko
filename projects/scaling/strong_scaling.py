import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

#import pymc3 as mcmc
import emcee
from scipy.optimize import minimize

from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import colorbar

import re
from cycler import cycler

def read_out_file(fname):
    data = {}
    data['laps'] = []

#--- add_cur                 0.062%  |  time:  2.44134 ms /100  (0.0024413443)
    #REcomp  = re.compile(r'---\s+(\S+)\s\S*(\S+)')
    REcomp  = re.compile(r'---\s+(\S+).*\((\S+)\)')
    RElap = re.compile(r'------ lap: (\S+)')

    lines = open(fname).readlines()
    for line in lines:
        line = line.strip()
        
        m = REcomp.match(line)
        if m:
            key = m.group(1)
            val = float(m.group(2))

            if key in data:
                data[key].append(val)
            else:
                data[key] = [val]

        m2 = RElap.match(line)
        if m2:
            data['laps'].append(float( m2.group(1) ))

    for k in data:
        data[k] = np.array(data[k][1:])

    return data



if __name__ == "__main__":
    fig = plt.figure(1, figsize=(3.487, 3.0))

    plt.rc('font',  family='serif')
    plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.3)
    #gs.update(wspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    
    for ax in axs:
        ax.minorticks_on()

    axs[0].set_yscale('log')
    axs[0].set_xscale('log')

    #axs[0].set_ylim((1.0e-9, 3.0e-6))
    #axs[0].set_xlim((5.0e-4, 1.0e+1))  #A
    axs[0].set_xlim((100, 1e4))
    axs[0].set_ylim((1e-3, 2e1))
    #axs[0].set_ylim((1.0e-9, 3.0e-6))

    axs[0].set_xlabel(r"Number of cores")
    axs[0].set_ylabel(r"Lap step time (s)")

    colcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    #colcycle = colcycle[0:5]

    custom_cycler = (
            cycler(linestyle=['solid', 'dashed', 'dotted','dashdot']) *
            #cycler(color=['skyblue', 'magenta', 'darkorange', 'y']) *
            cycler(color=colcycle) *
            cycler(lw=[1., 1.0]) 
                 )
    axs[0].set_prop_cycle(custom_cycler)



    filename = "sc3_c256_x512.out"
    data256 = read_out_file(filename)
    data256['procs'] = 256
    data256['nx'] = 5120**2
    data256['ppc'] = 64*2

    filename = "sc35_c512_x512.out"
    data512 = read_out_file(filename)
    data512['procs'] = 512
    data512['nx'] = 5120**2
    data512['ppc'] = 64*2

    #filename = "sc3b_c1024_x512.out"
    filename = "sc375_c1024_x512.out"
    data1024 = read_out_file(filename)
    data1024['procs'] = 1024
    data1024['nx'] = 5120**2
    data1024['ppc'] = 64*2

    filename = "sc39_c2048_x512.out"
    data2048 = read_out_file(filename)
    data2048['procs'] = 2048
    data2048['nx'] = 5120**2
    data2048['ppc'] = 64*2

    filename = "sc4_c4096_x512.out"
    data4096 = read_out_file(filename)
    data4096['procs'] = 4096
    data4096['nx'] = 5120**2
    data4096['ppc'] = 64*2


    #re-organize data as a function of procs per routine
    wdata = {}
    wdata['procs'] = [256, 512, 1024, 2048, 4096]
    wdata['total'] = []

    skip_keys = ['total', 'laps', 'io', 'init','step','avg:','std:','procs','nx','ppc','norm']
    for data in [data256, data512, data1024, data2048, data4096]:
        print(data['procs'])

        #data['norm'] = data['procs']/(data['nx']*data['ppc'])
        #data['norm'] = 1.0/(data['nx']*data['ppc'])
        data['norm'] = 1.0
        total = np.zeros(len(data['laps']))

        for key in data:
            if key in skip_keys:
                continue
            total[:] += data[key]*data['norm']

            #if np.max(data['norm']*data[key]) < 1.0e-8:
            #    continue
            #print(key)

            labelk = key.replace("_", "\_")
            labelk = "\\texttt{"+labelk+"}"

            proc = data['procs']

            avg = np.mean(data[key]*data['norm'])
            #axs[0].plot([proc], [avg], label=labelk)

            if key in wdata:
                wdata[key].append(avg)
            else:
                wdata[key] = [avg]

        tavg  = np.mean(total)
        #axs[0].plot([proc], [tavg], "k.", label=labelk)
        wdata['total'].append(tavg)


    skip_keys2 = ['procs']
    for key in wdata:
        if key in skip_keys2:
            continue

        #if np.max(wdata[key]) < 1.0e-8:
        #    continue
        print(key)
        print(wdata[key])

        labelk = key.replace("_", "\_")
        labelk = "\\texttt{"+labelk+"}"

        if key == 'total':
            bl, = axs[0].plot(wdata['procs'], wdata[key], "k", label=labelk)
        else:
            bl, = axs[0].plot(wdata['procs'], wdata[key], label=labelk)
        axs[0].plot(wdata['procs'], wdata[key], ".", color=bl.get_color() )


    cmin=128.
    cmax=16384.
    t0 = 17.
    t1 = t0*(cmin/cmax)
    axs[0].plot([cmin, cmax], [t0, t1], "k:", linewidth=2.5, alpha=0.6)


    handles, labels = axs[0].get_legend_handles_labels()

    axs[0].legend(handles, labels, fontsize=6)

    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), 
            loc="lower left",
            mode="expand", 
            borderaxespad=0, 
            ncol=3,
            fontsize=5,
            )

    #--------------------------------------------------
    axleft    = 0.18
    axbottom  = 0.15
    axright   = 0.96
    axtop     = 0.63
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)
    fname = 'strong_scaling.pdf'
    plt.savefig(fname)

