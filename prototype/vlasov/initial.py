import numpy as np

def initial(prm):

    #coordinates
    xx = np.arange(-3, prm.nx+3) * prm.dx
    vx = np.zeros((prm.nv+6, prm.ns))
    for kk in range(prm.ns):
        vx[:, kk] = np.linspace(prm.vmin[kk], prm.vmax[kk], prm.nv+6)
    #print xx,vx


    for kk in range(prm.ns):
        prm.dv[kk] = vx[1, kk] - vx[0, kk]
        prm.qn[kk] = prm.dv[kk]/( prm.qm[kk] * prm.wp[kk]**2 * prm.dx ) 

    kx = np.zeros(prm.nx)
    kv = np.zeros(prm.nv + 6)


    nkx = np.int( np.floor(prm.nx*0.5) )
    nkv = np.int( np.floor(prm.nv*0.5) + 3 )

    kx[0:nkx] = np.arange(0, nkx)/(prm.nx * prm.dx)*2.0*np.pi
    kv[0:nkv] = np.arange(0, nkv)/(prm.nv + 6)*2.0*np.pi
    for ii in range(nkx, prm.nx):
        kx[ii] = -kx[2*nkx + 2 - ii]

    for jj in range(nkv, prm.nv+6):
        kv[jj] = -kv[2*nkv +2 - jj]
    

    #field initialization
    fex = np.zeros(prm.nx + 6) #full integer grid
    ajx = np.zeros(prm.nx + 6)
    #ex  = (fex[0:prm.nx+5] + fex[1:prm.nx+6])/2.0 #half integer grid
    ex = np.zeros(prm.nx + 5)
    for ii in range(0, prm.nx+5):
        ex[ii] = (fex[ii] + fex[ii+1])/2.0


    #particle initialization
    ff = np.zeros( (prm.nv+6, prm.nx+6, prm.ns) )
    gx = np.zeros( (prm.nv+6, prm.nx+6, prm.ns) )
    gv = np.zeros( (prm.nv+6, prm.nx+6, prm.ns) )

    gam = 3
    wpe = np.sqrt( np.sum( prm.wp**2 * (-prm.qm) ) )

    for kk in range(prm.ns):
        if prm.qm[kk] < 0:
            ww = np.sqrt( wpe**2 + gam*prm.vt[kk]**2 * kx**2 )
        else:
            ww = prm.vd[kk]*kx
            prm.noise = 3

        #noise 
        prm.nmode = nkx
        amp = max( prm.pamp, prm.namp )
        prm.pamp = amp/prm.nx
        prm.namp = amp/prm.nx

        pphs = np.random.rand(prm.nmode)*360.0
        nphs = np.random.rand(prm.nmode)*360.0

        dn_noise = np.ones(prm.nx+6) #full integer grids
        dd_noise = np.zeros(prm.nx+6)
        vd_noise = np.ones(prm.nx+6) * prm.vd[kk]
        vt_noise = np.ones(prm.nx+6) * prm.vt[kk]

        for ll in range(prm.nmode):
            dn_noise = dn_noise - prm.pamp * np.sin(-kx[1 + ll]*xx + pphs[ll]/180*np.pi) * kx[1 + ll] \
                                + prm.namp * np.sin( kx[1 + ll]*xx + nphs[ll]/180*np.pi) * kx[1 + ll]
        
            dd_noise = dd_noise - prm.pamp * np.cos(-kx[1 + ll]*xx + pphs[ll]/180*np.pi) * kx[1 + ll]**2 \
                                + prm.namp * np.cos( kx[1 + ll]*xx + nphs[ll]/180*np.pi) * kx[1 + ll]**2


            vd_noise = vd_noise - prm.pamp * np.sin(-kx[1 + ll]*xx + pphs[ll]/180*np.pi)*(ww[ll] - prm.vd[kk]*kx[1 + ll]) \
                                - prm.namp * np.sin(-kx[1 + ll]*xx + nphs[ll]/180*np.pi)*(ww[ll] + prm.vd[kk]*kx[1 + ll])


        for ii in range(prm.nx+6):
            for jj in range(prm.nv+6):
                ff[jj, ii, kk] = np.exp(-(vx[jj, kk] - vd_noise[ii])**2/(2*vt_noise[ii]**2)) \
                / (np.sqrt(2*np.pi)*vt_noise[ii])*dn_noise[ii]

                gx[jj, ii, kk] = np.exp(-(vx[jj, kk] - vd_noise[ii])**2/(2*vt_noise[ii]**2)) \
                / (np.sqrt(2*np.pi)*vt_noise[ii])*dn_noise[ii] * prm.dx
                
                gv[jj, ii, kk] = -np.exp(-(vx[jj, kk] - vd_noise[ii])**2/(2*vt_noise[ii]**2)) \
                / (np.sqrt(2*np.pi)*vt_noise[ii])*dn_noise[ii] \
                * (vx[jj, kk] - vd_noise[ii])/(vt_noise[ii]**2) * prm.dv[kk]


    return ff, gx, gv, ex, fex, ajx, xx, vx, kx, kv





