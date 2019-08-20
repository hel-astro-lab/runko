import numpy as np
import matplotlib.pyplot as plt
from combine_files import get_file_list
from parser import parse_input

from plot_shock import TestParticles

if __name__ == "__main__":
    #Create Container and direct to file
    prtcls = TestParticles()
    conf, fdir, args = parse_input()

    fdir += '/'
    fname_F = "flds"
    fname_A = "analysis"
    fname_P = "test-prtcls"

    #Add all files to container
    files_F = get_file_list(fdir, fname_F)
    for lap, f in enumerate(files_F):
        info = {}
        info['lap'] = lap*conf.interval
        info['fields_file'  ]   = files_F[lap]
        info['skindepth'] = conf.c_omp/conf.stride
    
        info['particle_file'] = fdir + 'test-prtcls'
        prtcls.add_file(info)
    
    #Read file with top 10 particles to obtain keys
    f = open('{}/10_prtcls.txt'.format(conf.outdir), "r")
    contents = f.readlines()
    topkeys_id = []
    topkeys_proc = []
    topkeys = []

    j = 0
    for i in contents:
        topkeys_id.append(i.split(',')[0])
        topkeys_proc.append(i.split(',')[1])
        topkeys.append((int(topkeys_id[j]),int(topkeys_proc[j])))
        j += 1
    f.close()
    print(topkeys)
    print("Keys loaded...")

    plt.fig = plt.figure(1, figsize=(3.54, 3), dpi=300)
    plt.rcParams.update({'font.size': 8})
    plt.subplots_adjust(wspace = 0, top = 0.98, bottom = 0.15, left = 0.2, right = 0.95)
    gs = plt.GridSpec(1, 2)
    ax = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])

    ax.set_xlabel('$x$ ($c/\omega_p$)', color = 'black')
    ax.set_ylabel('$t$ ($\omega_{p}^{-1}$)', color = 'black')
    ax.tick_params(axis='x', colors='black')
    ax.tick_params(axis='y', colors='black')
    ax.minorticks_on()
    
    ax2.set_xlabel(r'$\gamma$', color = 'black')
    ax2.tick_params(axis='x', colors='black')
    ax2.minorticks_on()
    ax2.set_yticklabels([])
    ax2.set_xscale('log')

    Nolaps = int((conf.Nt/conf.interval)+1)
    laplen = range(Nolaps)
    laplen = np.array(laplen)
    t_omp = laplen*(conf.cfl/conf.c_omp)
    
    c = 0
    for key in topkeys:
        p = prtcls.prtcls[key]
        p_np = np.array(p)
        x = p_np[:,0]
        
        ax.plot(x, t_omp)
        
        ux = p_np[:,3]
        uy = p_np[:,4]
        uz = p_np[:,5]
    
        gamma = (1 + ux**2 + uy**2 + uz**2)**0.5
        ax2.plot(gamma, laplen)
        c += 1
        if c == 5:
            break
    
    #plt.tight_layout()
    plt.savefig('{}/prtcl_path.pdf'.format(conf.outdir))
    #omp = conf.cfl/conf.c_omp #Plasma frequency
    #qe = (omp**2.*conf.gamma)/((conf.ppc*.5)*(1.+abs(conf.me/conf.mi)))
    #print((conf.ppc * abs(qe) / (conf.gamma ))**0.5)
    
    

