import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

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
        data[k] = np.array(data[k][0:])

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

    #axs[0].set_yscale('log')
    #axs[0].set_xscale('log')

    #axs[0].set_ylim((1.0e-9, 3.0e-6))
    #axs[0].set_xlim((5.0e-4, 1.0e+1))  #A

    axs[0].set_ylim((0, 2.0))
    axs[0].set_xlim((0, 40))

    axs[0].set_xlabel(r"Tile size ($N_x^3$)")
    axs[0].set_ylabel(r"Absolute time per lap (s)")

    colcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    #colcycle = colcycle[0:5]

    custom_cycler = (
            cycler(linestyle=['solid', 'dashed', 'dotted','dashdot']) *
            #cycler(color=['skyblue', 'magenta', 'darkorange', 'y']) *
            cycler(color=colcycle) *
            cycler(lw=[1., 1.0]) 
                 )
    axs[0].set_prop_cycle(custom_cycler)

    filename = "8.out"
    data0 = read_out_file(filename)
    data0['tile'] = 8

    filename = "15.out"
    data1 = read_out_file(filename)
    data1['tile'] = 15

    filename = "30.out"
    data2 = read_out_file(filename)
    data2['tile'] = 30


    #re-organize data as a function of procs per routine
    wdata = {}
    wdata['tile'] = [8, 15, 30]
    wdata['total'] = []

    skip_keys = ['total', 'laps', 'io', 'init','step','avg:','std:','procs','nx','ppc','norm', 'tile', 'mpi_e0']
            #'sliding_box', 'rwall', 'walls',]

    for data in [data0, data1, data2, ]:
        print(data['tile'])
        print(data)

        lent = len(data['sliding_box']) # get number of laps processed
        total = np.zeros(lent)

        for key in data:
            if key in skip_keys:
                continue

            labelk = key.replace("_", "\_")
            labelk = "\\texttt{"+labelk+"}"

            if False: # calculate mean from given times
                total[:] += data[key][:lent]
                avg = np.mean(data[key][:lent])
                #axs[0].plot([proc], [avg], label=labelk)

                if key in wdata:
                    wdata[key].append(avg)
                else:
                    wdata[key] = [avg]


            if True: # pick one time slice; here first lap
                ind = 0
                total[:] += data[key][ind]
                avg = data[key][ind]

                #axs[0].plot([proc], [avg], label=labelk)

                if key in wdata:
                    wdata[key].append(avg)
                else:
                    wdata[key] = [avg]


        tavg  = np.mean(total)
        #axs[0].plot([proc], [tavg], "k.", label=labelk)
        wdata['total'].append(tavg)


    print(wdata['total'])

    skip_keys2 = ['procs', 'tile']
    for key in wdata:
        if key in skip_keys2:
            continue

        #if np.max(wdata[key]) < 1.0e-8:
        #    continue
        #print(key)
        #print(wdata[key])

        labelk = key.replace("_", "\_")
        labelk = "\\texttt{"+labelk+"}"

        if key == 'total':
            bl, = axs[0].plot(wdata['tile'], wdata[key], "k", label=labelk)
        else:
            bl, = axs[0].plot(wdata['tile'], wdata[key], label=labelk)
        axs[0].plot(wdata['tile'], wdata[key], ".", color=bl.get_color() )



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
    fname = 'tile_scaling.pdf'
    plt.savefig(fname)

