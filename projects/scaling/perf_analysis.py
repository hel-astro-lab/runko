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
    new_lap_mode = False

    data['laps'] = []

#--- add_cur                 0.062%  |  time:  2.44134 ms /100  (0.0024413443)
    #REcomp  = re.compile(r'---\s+(\S+)\s\S*(\S+)')
    REcomp  = re.compile(r'---\s+(\S+).*\((\S+)\)')
    RElap = re.compile(r'------ lap: (\S+)')

    lines = open(fname).readlines()
    for line in lines:
        line = line.strip()

        #if line == "--------------------------------------------------":
        #    #activate if not on
        #    if new_lap_mode:
        #        new_lap_mode = False
        #    else:
        #        new_lap_mode = True
        #if new_lap_mode:
        #    print(line)
        
        m = REcomp.match(line)
        if m:
            key = m.group(1)
            val = float(m.group(2))

            if key in data:
                data[key].append(val)
            else:
                data[key] = []

        m2 = RElap.match(line)
        if m2:
            data['laps'].append(float( m2.group(1) ))

    data['laps'].pop()

    for k in data:
        data[k] = np.array(data[k])

    return data



if __name__ == "__main__":
    fig = plt.figure(1, figsize=(3.487, 2.5))

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
    #axs[0].set_xscale('log')

    axs[0].set_ylim((1.0e-9, 3.0e-6))
    #axs[0].set_xlim((5.0e-4, 1.0e+1))  #A

    axs[0].set_xlabel(r"laps $n$")
    axs[0].set_ylabel(r"Push time (s prctl$^{-1}$ core$^{-1}$)")

    colcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    #colcycle = colcycle[0:5]

    custom_cycler = (
            cycler(linestyle=['solid', 'dashed', 'dotted','dashdot']) *
            #cycler(color=['skyblue', 'magenta', 'darkorange', 'y']) *
            cycler(color=colcycle) *
            cycler(lw=[1.]) 
                 )
    axs[0].set_prop_cycle(custom_cycler)

    if False:
        filename = "64_16.out"
        data = read_out_file(filename)

        data['procs'] = 16
        data['nx'] = 640**2
        data['ppc'] = 64*2

    if True:
        filename = "128_64.out"
        data = read_out_file(filename)

        data['procs'] = 128
        data['nx'] = 1280**2
        data['ppc'] = 64*2

    if False:
        filename = "256_256.out"
        data = read_out_file(filename)

        data['procs'] = 256
        data['nx'] = 2560**2
        data['ppc'] = 64*2

    if False:
        filename = "512_1024.out"
        data = read_out_file(filename)

        data['procs'] = 1024
        data['nx'] = 5120**2
        data['ppc'] = 64*2

    data['norm'] = data['procs']/(data['nx']*data['ppc'])

    total = np.zeros(len(data['laps']))

    skip_keys = ['total', 'laps', 'io', 'init','step','avg:','std:','procs','nx','ppc','norm']

    #first get total
    for key in data:
        if key in skip_keys:
            continue
        total[:] += data[key]*data['norm']
    axs[0].plot(data['laps'], total, "k-", label="total")

    #then component-wise
    for key in data:
        if key in skip_keys:
            continue

        #print(data[key])

        #if np.max(data['norm']*data[key]) < 1.0e-8:
        #    continue

        print("plotting ", key)
        labelk = key.replace("_", "\_")
        labelk = "\\texttt{"+labelk+"}"
        axs[0].plot(data['laps'], data['norm']*data[key], label=labelk)

    print("mean total push time: ", np.mean(total))
    print("min total push time: ", np.min(total))
    print("max total push time: ", np.max(total))

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
    fname = 'perf_analysis.pdf'
    plt.savefig(fname)

