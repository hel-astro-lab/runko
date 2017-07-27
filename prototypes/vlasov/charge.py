import numpy as np

#charge density
def charge(ff, prm):
    rhos = np.zeros( (prm.nx + 6, prm.ns) )
    
    #sum over velocity dimension
    for kk in range(prm.ns):
        for ii in range(prm.nx + 6):
            rhos[ii, kk] = np.sum( ff[:, ii, kk] ) * prm.qn[kk]

    #sum over spatial dimension
    rho = np.sum(rhos, 1)

    # XXX nx+3 or nx+2?
    rho0 = np.mean( rho[3:prm.nx+3] )
    rho = rho - rho0

    return rho
