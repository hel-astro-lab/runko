import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import h5py as h5
import sys, os

import argparse

from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import colorbar

import copy


if __name__ == "__main__":
    fig = plt.figure(1, figsize=(3.487, 10.0))
    #fig = plt.figure(1, figsize=(7.74, 10.0))
    

    plt.rc('font',  family='sans-serif')
    #plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.0)
    gs.update(wspace = 0.0)
    
    #XXX add panels here
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    #axs.append( plt.subplot(gs[1,0]) )

    
    for ax in axs:
        ax.minorticks_on()

    if len(axs) > 1:
        axs[0].tick_params(which='x', direction="in")
        axs[0].axes.get_xaxis().set_ticks([])

    #--------------------------------------------------
    if False:
        left, bottom, width, height = [0.65, 0.67, 0.3, 0.15]
        ax2 = fig.add_axes([left, bottom, width, height])
        #ax2.set_yscale('log')
        #ax2.set_xscale('log')
        ax2.set_xlabel(r"$t$ $(l_0/c)$", fontsize=5)
        ax2.set_ylabel(r"$\theta$", fontsize=5)
        ax2.xaxis.set_tick_params(labelsize=4)
        ax2.yaxis.set_tick_params(labelsize=4)
        ax2.tick_params(which='both', direction="in")


    #--------------------------------------------------
    # set up plot
    tmin = -0.6
    tmax =  1.3

    ymin = -1
    global ycur
    ycur = 0
    ydt = 1.0

    #ax.set_yscale('log')
    #ax.set_xscale('log')

    ax.set_xlim((tmin, tmax))
    #ax.set_xlim((-1.5, 1.5))
                           
    ax.set_xlabel(r'Time $t^n$ $(\Delta t)$')
    ax.set_ylabel(r"Sub-cycles")


    #variables = \
    #{
    #        'J':{'in':-0.5, 'in_bc':-0.5, 'out_bc':-0.5, 'vir':-0.5, 'ext':-0.5},
    #        'B':{'in':-0.0, 'in_bc':-0.0, 'out_bc':-0.0, 'vir':-0.0, 'ext':-0.0},
    #        'E':{'in':-0.5, 'in_bc':-0.5, 'out_bc':-0.5, 'vir':-0.5, 'ext':-0.5},
    #        'x':{'in':-0.0, 'in_bc':-0.0, 'out_bc':-0.0, 'vir':-0.0, 'ext':-0.0},
    #        'v':{'in':-0.5, 'in_bc':-0.5, 'out_bc':-0.5, 'vir':-0.5, 'ext':-0.5},
    #}

    #opposite staggering
    variables = \
    {
            'J':{'in':-0.0, 'in_bc':-0.0, 'out_bc':-0.0, 'vir':-0.0, 'ext':-0.0},
            'B':{'in':-0.5, 'in_bc':-0.5, 'out_bc':-0.5, 'vir':-0.5, 'ext':-0.5},
            'E':{'in':-0.0, 'in_bc':-0.0, 'out_bc':-0.0, 'vir':-0.0, 'ext':-0.0},
            'x':{'in':-0.0, 'in_bc':-0.0, 'out_bc':-0.0, 'vir':-0.0, 'ext':-0.0},
            'v':{'in':-0.5, 'in_bc':-0.5, 'out_bc':-0.5, 'vir':-0.5, 'ext':-0.5},
    }

    def get_tile_type(regime):
        ts = []
        if regime == 'local':
            ts.append('in')
            ts.append('ext') 
        if regime == 'all':
            ts.append('in')
            ts.append('vir')
            ts.append('ext') 
        if regime == 'vir': 
            ts.append('vir')
        if regime == 'mpi': 
            ts.append('vir')
        if regime == 'bc':
            ts.append('in_bc')
            ts.append('out_bc')

        return ts
    

    #visualize_vars = ['in', 'vir','in_bc','out_bc','ext']
    visualize_vars = ['in', 'vir','in_bc','out_bc']
    #visualize_vars = ['in']


    #shift variable in time completely; done to counter sub-cycle shifts
    def shift_variable(var, dt):
        variables[var]['in']    += dt
        variables[var]['in_bc'] += dt
        variables[var]['out_bc']+= dt
        variables[var]['vir']   += dt
        variables[var]['ext']   += dt



    def apply_routine(algo):
        list_vars = algo['var']
        dt      = algo['dt']
        regime  = algo['regime']

        for var in list_vars:

            #get type of tiles/regions affected
            tile_types = get_tile_type(regime)
            for tt in tile_types:

                #plot initial condition
                if tt in visualize_vars:
                    t0,y0 = plot_var(ax, var, tt)

                #track changes
                torig = copy.deepcopy(variables[var][tt])

                if not(algo['copy']):
                    variables[var][tt] += dt
                else:
                    #special copies branch
                    
                    if regime == 'mpi':
                        print("copy mpi {}: {} <- {}".format(
                            var, 
                            variables[var]['vir'],
                            variables[var]['ext'],
                            ))

                        variables[var]['vir'] = copy.deepcopy(variables[var]['ext'])

                    if regime == 'bc':
                        if tt == 'in_bc':
                            #local tiles copy neighboring local tile time step
                            variables[var]['in_bc']  = copy.deepcopy(variables[var]['in'])
                        if tt == 'out_bc':
                            #outer tiles copy virtual tile time step
                            variables[var]['out_bc'] = copy.deepcopy(variables[var]['vir'])

                #track changes
                tend = copy.deepcopy(variables[var][tt])
                if not(torig == tend):

                    #plot final configuration
                    if tt in visualize_vars:
                        t1,y1 = plot_var(ax, var, tt)

                        #arrow (t0,y0) -> (t1,ycur-1)
                        padd = 0.03
                        arrowprops=dict(facecolor='black', width=0.5, headwidth=3.0, headlength=5.0, linewidth=0.2) #, shrink=1.0), )
                        plt.annotate('', xy=(tend-padd,ycur-ydt),  xytext=(torig+padd,y0),  arrowprops=arrowprops)

        ax.axhline(ycur-ydt/2, linestyle='dotted', alpha=0.4)

        #latex version
        if False:
            strname = algo['name'].replace('_','\_')
            #strname = r'${{{0}}}$'.format(strname)
            strname = "\\texttt{"+strname+"}"
        else:
            strname = algo['name']
        axs[0].text(tmax + 0.03, ycur-ydt, strname, va='center', ha='left', fontsize=6, fontfamily='monospace')
    


    def plot_var(ax, var, tt, pos=None):
        global ycur
        
        def get_style(tt):
            style = {}

            #bbox=dict(facecolor='red', alpha=0.5, boxstyle='square'))
            if tt == 'in':
                #style = {'color':'black'}
                style = {
                        'color':'black', 
                        'bbox':{'facecolor':'gray', 'alpha':0.1, 'boxstyle':'square', 'pad':0.12},
                        }
            if tt == 'vir':
                style = {
                        'color':'red', 
                        'bbox':{'facecolor':'red', 'alpha':0.5, 'boxstyle':'circle', 'pad':0.08},
                        }
            if tt == 'in_bc':
                style = {
                        'color':'red', 
                        'bbox':{'facecolor':'blue', 'alpha':0.5, 'boxstyle':'square', 'pad':0.08},
                        }
            if tt == 'out_bc':
                style = {
                        'color':'red', 
                        'bbox':{'facecolor':'red', 'alpha':0.5, 'boxstyle':'square', 'pad':0.08},
                        }
            if tt == 'ext':
                style = {
                        'color':'green', 
                        #'bbox':{'facecolor':'red', 'alpha':0.5, 'boxstyle':'circle'},
                        }

            return style
    
        t0 = variables[var][tt]
        ss = get_style(tt)

        strvar = r'${{{}}}$'.format(var)
        ax.text(t0, ycur, strvar, ha='center', va='center',fontsize=7, **ss)

        y0 = copy.deepcopy(ycur)
        ycur += ydt

        return t0, y0

            
    #PIC loop
    if True:

        #--------------------------------------------------
        for lap in range(0,1):


            #update half b
            step1 = [
                    #{'name':'mpi_e1','var':['E', ]    , 'regime':'mpi',  'dt':0,   'copy':True  , },
                    #{'name':'upd_bc','var':['E',] ,     'regime':'bc',   'dt':0,   'copy':True , } ,
                    #{'name':'f_bc'  ,'var':['B','E']  , 'regime':'bc',   'dt':0,   'copy':False , }, # print loc 
                    #{'name':'f_bc'  ,'var':['B','E']  , 'regime':'local','dt':0,   'copy':False , }, # print loc
                    {'name':'push_b','var':['B',]     , 'regime':'local','dt':0.5, 'copy':False , },
                    {'name':'mpi_b1','var':['B',]     , 'regime':'mpi',  'dt':0,   'copy':True  , },
                    {'name':'upd_bc','var':['B',] ,     'regime':'bc',   'dt':0,   'copy':True , } ,
            ]
            for algo in step1:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)


            ##push prtcls
            step2 = [
                    #{'name':'f_bc'  ,    'var':['B','E'] , 'regime':'bc',    'dt':0, 'copy':False , }, #print 
                    #{'name':'f_bc'  ,    'var':['B','E']  ,'regime':'local', 'dt':0, 'copy':False , }, #print 
                    {'name':'interp_em', 'var':['E','B',], 'regime':'local', 'dt':0, 'copy':False},
                    {'name':'push',      'var':['v','x',], 'regime':'local', 'dt':1, 'copy':False},
                    ]
            for algo in step2:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)


            # push second half
            step3 = [
                    #{'name':'f_bc'  ,'var':['B','E'] , 'regime':'bc',    'dt':0, 'copy':False , }, #print 
                    #{'name':'f_bc'  ,'var':['B','E']  ,'regime':'local', 'dt':0, 'copy':False , }, #print 
                    {'name':'push_b','var':['B',],     'regime':'local','dt':0.5, 'copy':False,},
                    {'name':'mpi_b2','var':['B',]     , 'regime':'mpi',  'dt':0,   'copy':True  , },
                    {'name':'upd_bc','var':['B',] ,     'regime':'bc',   'dt':0,   'copy':True , } ,
                    ]
            for algo in step3:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)


            ##push e
            step4 = [
                    #{'name':'f_bc'  ,  'var':['B','E'] , 'regime':'bc',    'dt':0, 'copy':False , }           ,
                    #{'name':'f_bc'  ,  'var':['B','E'] , 'regime':'local', 'dt':0, 'copy':False , }           ,
                    {'name':'push_e',  'var':['E',],      'regime':'local','dt':1.0, 'copy':False,},
                    ]
            for algo in step4:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)


            #deposit
            step5 = [
                    #{'name':'j_bc'  ,    'var':['v'] ,     'regime':'bc',    'dt':0  , 'copy':False , }           ,
                    #{'name':'j_bc'  ,    'var':['v','x'] , 'regime':'local', 'dt':0  , 'copy':False , }           ,
                    {'name':'comp_curr', 'var':['J',], 'regime':'local', 'dt':1  , 'copy':False ,},
                    {'name':'mpi_curr',  'var':['J',], 'regime':'mpi',   'dt':0  , 'copy':True ,},
                    {'name':'exc_cur',   'var':['J'],  'regime':'bc',    'dt':0  , 'copy':True ,},
                    ]
            for algo in step5:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)


            ##filter
            step7 = [
                    {'name':'filter',  'var':['J',],      'regime':'local', 'dt':0.10  ,'copy':False ,},
                    {'name':'mpi_j0',  'var':['J',],      'regime':'mpi',   'dt':0  ,'copy':True ,},
                    {'name':'upd_bc0', 'var':['J',],      'regime':'bc',    'dt':0  ,'copy':True ,},
                    #{'name':'j_bc'  ,  'var':['J'],       'regime':'bc',    'dt':0,  'copy':False , }           ,
                    ]
            for algo in step7:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ##second pass
            #for algo in step7:
            #    apply_routine(algo)

            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ##shift back
            #shift_variable('J', -0.2)

            ##add current
            step8 = [
                    #{'name':'j_bc'  ,  'var':['J','E',],  'regime':'bc',    'dt':0, 'copy':False , }           ,
                    #{'name':'j_bc'  ,  'var':['J',],      'regime':'local', 'dt':0, 'copy':False , }           ,
                    {'name':'add_cur', 'var':['E',],      'regime':'local', 'dt':0.1,'copy':False ,},
                    {'name':'mpi_e2',  'var':['E',],      'regime':'mpi',   'dt':0  ,'copy':True ,},
                    {'name':'upd_bc',  'var':['E',] ,     'regime':'bc',    'dt':0,   'copy':True , } ,
                    ]
            for algo in step8:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            #shift back
            #shift_variable('E', -0.1)


            ##mpi prtcls
            step6 = [
                    {'name':'outg_prtcls',     'var':['x','v',], 'regime':'bc',   'dt':0  , 'copy':True ,},
                    {'name':'mpi_prtcls',      'var':['x','v'],   'regime':'mpi',  'dt':0  , 'copy':True ,},
                    {'name':'get_inc_prtcls',  'var':['x','v',],  'regime':'bc',   'dt':0  , 'copy':True ,},
                    {'name':'del_trnsf_prtcls','var':['x','v',], 'regime':'local',   'dt':0  , 'copy':False ,},
                    {'name':'del_vir_prtcls',  'var':['x','v',], 'regime':'vir',     'dt':0  , 'copy':False ,},
                    ]
            for algo in step6:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)






            #--------------------------------------------------

            #update half b
            #step1 = [
            #        {'name':'mpi_b1','var':['B', ]    , 'regime':'mpi',  'dt':0,   'copy':True  , },
            #        {'name':'upd_bc','var':['E','B',] , 'regime':'bc',   'dt':0,   'copy':True , } ,
            #        {'name':'f_bc'  ,'var':['B','E']  , 'regime':'bc',   'dt':0.0, 'copy':False , },
            #        {'name':'f_bc'  ,'var':['B','E']  , 'regime':'local',   'dt':0.0, 'copy':False , },
            #        {'name':'push_b','var':['B',]     , 'regime':'local','dt':0.5, 'copy':False , },
            #        {'name':'mpi_b2','var':['B',]     , 'regime':'mpi',  'dt':0,   'copy':True  , },
            #        {'name':'upd_bc','var':['E','B',] , 'regime':'bc',   'dt':0,   'copy':True , } ,
            #]
            #for algo in step1:
            #    apply_routine(algo)
            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ##push prtcls
            #step2 = [
            #        {'name':'f_bc'  ,    'var':['B','E'] , 'regime':'bc',    'dt':0, 'copy':False , }           ,
            #        {'name':'interp_em', 'var':['E','B',], 'regime':'local', 'dt':0, 'copy':False},
            #        {'name':'push',      'var':['v','x',], 'regime':'local', 'dt':1, 'copy':False},
            #        ]
            #for algo in step2:
            #    apply_routine(algo)
            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            #step3 = [
            #        {'name':'push_b', 'var':['B',],     'regime':'local','dt':0.5, 'copy':False,},
            #        ]
            #for algo in step3:
            #    apply_routine(algo)
            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ##push e
            #step4 = [
            #        {'name':'mpi_e1',  'var':['B',],     'regime':'mpi',  'dt':0  , 'copy':True ,},
            #        {'name':'upd_bc2', 'var':['E','B',], 'regime':'bc',   'dt':0  , 'copy':True ,},
            #        {'name':'f_bc'  ,  'var':['B','E'] , 'regime':'bc',    'dt':0, 'copy':False , }           ,
            #        {'name':'f_bc'  ,  'var':['B','E'] , 'regime':'local', 'dt':0, 'copy':False , }           ,
            #        {'name':'push_e',  'var':['E',],      'regime':'local','dt':1.0, 'copy':False,},
            #        ]
            #for algo in step4:
            #    apply_routine(algo)
            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ##current
            #step5 = [
            #        {'name':'j_bc'  ,    'var':['v'] , 'regime':'bc',    'dt':0  , 'copy':False , }           ,
            #        {'name':'j_bc'  ,    'var':['v','x'] , 'regime':'local', 'dt':0  , 'copy':False , }           ,
            #        {'name':'comp_curr', 'var':['J',], 'regime':'local', 'dt':1  , 'copy':False ,},
            #        {'name':'mpi_curr',  'var':['J',], 'regime':'mpi',   'dt':0  , 'copy':True ,},
            #        {'name':'exc_cur',   'var':['J'],  'regime':'bc',    'dt':0  , 'copy':True ,},
            #        ]
            #for algo in step5:
            #    apply_routine(algo)
            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)


            ##mpi prtcls
            #step6 = [
            #        {'name':'outg_prtcls', 'var':['x','v',], 'regime':'bc',   'dt':0  , 'copy':True ,},
            #        {'name':'mpi_prtcls', 'var':['x','v'],   'regime':'mpi',  'dt':0  , 'copy':True ,},
            #        {'name':'get_inc_prtcls', 'var':['x','v',],  'regime':'bc',   'dt':0  , 'copy':True ,},
            #        {'name':'del_trnsf_prtcls', 'var':['x','v',], 'regime':'local',   'dt':0  , 'copy':False ,},
            #        {'name':'del_vir_prtcls',   'var':['x','v',], 'regime':'vir',     'dt':0  , 'copy':False ,},
            #        ]
            #for algo in step6:
            #    apply_routine(algo)
            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ##filter
            ##add current
            #step7 = [
            #        {'name':'mpi_j0',  'var':['J',],      'regime':'mpi',   'dt':0  ,'copy':True ,},
            #        {'name':'upd_bc0', 'var':['J',],      'regime':'bc',   'dt':0  , 'copy':True ,},
            #        {'name':'j_bc'  ,  'var':['J'],       'regime':'bc',    'dt':0,  'copy':False , }           ,
            #        {'name':'filter',  'var':['J',],      'regime':'local', 'dt':0.10  ,'copy':False ,},
            #        ]
            #for algo in step7:
            #    apply_routine(algo)
            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ##second pass
            #for algo in step7:
            #    apply_routine(algo)

            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ##shift back
            #shift_variable('J', -0.2)



            ##add current
            #step8 = [
            #        {'name':'j_bc'  ,  'var':['J','E',],  'regime':'bc',    'dt':0, 'copy':False , }           ,
            #        {'name':'j_bc'  ,  'var':['J',],      'regime':'local', 'dt':0, 'copy':False , }           ,
            #        {'name':'add_cur', 'var':['E',],      'regime':'local', 'dt':0.1,'copy':False ,},
            #        {'name':'mpi_e2',  'var':['E',],      'regime':'mpi',   'dt':0  ,'copy':True ,},
            #        ]
            #for algo in step8:
            #    apply_routine(algo)
            #ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ##shift back
            #shift_variable('E', -0.1)


            ##print out of current state
            #ax.axhline(ycur-ydt/2, linestyle='solid', alpha=0.8)

            #step9 = [
            #        {'name':'io_halos', 'var':['B', 'E', 'J', 'x','v',], 'regime':'bc',    'dt':0  , 'copy':False ,},
            #        {'name':'io_virs',  'var':['B', 'E', 'J', 'x','v',], 'regime':'vir',   'dt':0  , 'copy':False ,},
            #        {'name':'io_tiles', 'var':['B', 'E', 'J', 'x','v',], 'regime':'local', 'dt':0  , 'copy':False ,},
            #        ]
            #for algo in step9:
            #    apply_routine(algo)

            #ax.axhline(ycur-ydt/2, linestyle='solid', alpha=0.8)


    #FFE loop
    else:
        variables = \
        {
                'J' :{'in':-0.0, 'in_bc':-0.0, 'out_bc':-0.0, 'vir':-0.0, 'ext':-0.0},
                'J1':{'in':+0.0, 'in_bc':+0.0, 'out_bc':+0.0, 'vir':+0.0, 'ext':+0.0},
                'B' :{'in':-0.5, 'in_bc':-0.5, 'out_bc':-0.5, 'vir':-0.5, 'ext':-0.5},
                'E' :{'in':-0.0, 'in_bc':-0.0, 'out_bc':-0.0, 'vir':-0.0, 'ext':-0.0},
        }

        #--------------------------------------------------
        for lap in range(0,1):

            #push b
            step1 = [
                    #{'name':'f_bc'  ,'var':['B','E']  , 'regime':'bc',   'dt':0.0, 'copy':False , },
                    #{'name':'f_bc'  ,'var':['B','E']  , 'regime':'local','dt':0.0, 'copy':False , },
                    {'name':'push_b','var':['B',]     , 'regime':'all','dt':0.5, 'copy':False , },
            ]
            for algo in step1:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            #drift perp current
            step3 = [
                    #{'name':'mpi_b2','var':['B',]        , 'regime':'mpi',  'dt':0,   'copy':True  , },
                    {'name':'upd_bc','var':['E','B',]    , 'regime':'bc',   'dt':0,   'copy':True , } ,
                    {'name':'f_bc'  ,   'var':['B','E']  , 'regime':'bc',   'dt':0.0, 'copy':False , },
                    #{'name':'f_bc'  ,   'var':['B','E']  , 'regime':'local','dt':0.0, 'copy':False , },
                    {'name':'drift_cur','var':['J',]     , 'regime':'all','dt':0.1, 'copy':False , },
                    {'name':'deposit_j','var':['E',]     , 'regime':'all','dt':0.1, 'copy':False , },
            ]
            for algo in step3:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)
            shift_variable('J', -0.1)

            # parallel current
            step4 = [
                    {'name':'upd_bc','var':['E','B',]    , 'regime':'bc',   'dt':0,   'copy':True , } ,
                    {'name':'f_bc'  ,   'var':['B','E']  , 'regime':'bc',   'dt':0.0, 'copy':False , },
                    #{'name':'f_bc'  ,   'var':['B','E']  , 'regime':'local','dt':0.0, 'copy':False , },
                    {'name':'para_cur', 'var':['J1',]    , 'regime':'all','dt':0.1, 'copy':False , },
                    {'name':'deposit_j','var':['E',]     , 'regime':'all','dt':0.1, 'copy':False , },
            ]
            for algo in step4:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)
            shift_variable('J1', -0.1)

            #limiter
            step4 = [
                    {'name':'upd_bc','var':['E','B',]    , 'regime':'bc',   'dt':0,   'copy':True , } ,
                    {'name':'f_bc'  ,   'var':['B','E']  , 'regime':'bc',   'dt':0.0, 'copy':False , },
                    #{'name':'f_bc'  ,   'var':['B','E']  , 'regime':'local','dt':0.0, 'copy':False , },
                    {'name':'limiter',  'var':['E']  , 'regime':'local','dt':0.1, 'copy':False , },
            ]
            for algo in step4:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            shift_variable('E', -0.3)


            #push b
            step1 = [
                    {'name':'mpi_e1','var':['E',]        , 'regime':'mpi',  'dt':0,   'copy':True  , },
                    {'name':'upd_bc','var':['E','B',]    , 'regime':'bc',   'dt':0,   'copy':True , } ,
                    #{'name':'f_bc'  ,'var':['B','E']  , 'regime':'bc',   'dt':0.0, 'copy':False , },
                    #{'name':'f_bc'  ,'var':['B','E']  , 'regime':'local','dt':0.0, 'copy':False , },
                    {'name':'push_b','var':['B',]     , 'regime':'local',   'dt':0.5, 'copy':False , },
            ]
            for algo in step1:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            #push e
            step2 = [
                    {'name':'mpi_b1','var':['B',]        , 'regime':'mpi',  'dt':0,   'copy':True  , },
                    {'name':'upd_bc','var':['E','B',]    , 'regime':'bc',   'dt':0,   'copy':True , } ,
                    #{'name':'f_bc'  ,'var':['B','E']  , 'regime':'bc',   'dt':0.0, 'copy':False , },
                    #{'name':'f_bc'  ,'var':['B',]     , 'regime':'local','dt':0.0, 'copy':False , },
                    {'name':'push_e','var':['E',]     , 'regime':'all','dt':1.0, 'copy':False , },
                    #{'name':'mpi_e1','var':['E', ]    , 'regime':'mpi',  'dt':0,   'copy':True  , },
                    {'name':'upd_bc','var':['E','B',] , 'regime':'bc',   'dt':0,   'copy':True , } ,
            ]
            for algo in step2:
                apply_routine(algo)
            ax.axhline(ycur-ydt/2, linestyle='dashed', alpha=0.6)

            ax.axhline(ycur-ydt/2, linestyle='solid', alpha=0.8)


            #step drift current one time step (artificial since they are computed in place)
            shift_variable('J',  1.0)
            shift_variable('J1', 1.0)


    ax.set_ylim((ymin, ycur))

    #--------------------------------------------------    
    axleft    = 0.08
    axbottom  = 0.06
    axright   = 0.82
    axtop     = 0.98

    pos1 = axs[0].get_position()
    axwidth  = axright - axleft
    axheight = (axtop - axbottom)*0.02
    axpad = 0.01

    #cax = fig.add_axes([axleft, axtop + axpad, axwidth, axheight])
    #cb1 = colorbar.ColorbarBase(
    #        cax,
    #        cmap=cmap,
    #        norm=norm,
    #        orientation='horizontal',
    #        ticklocation='top')
    #cb1.set_label(r'$t$ $(l_0/c)$')

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    fdir = ''
    fname = fdir+'timestepping.pdf'
    plt.savefig(fname)

