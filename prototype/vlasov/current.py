import numpy as np

def current(ff, vx, prm):

    ajxs = np.zeros( (prm.nx+6, prm.ns) )
    for kk in range(prm.ns):
        aa = vx[:, kk] / prm.dx * prm.dt

        for ii in range(prm.nx+6):
            ajxs[ii, kk] = np.sum( ff[:, ii, kk] ) * prm.qn[kk]

    ajx = np.sum(ajxs, 1)
    ajx -= np.mean( ajx[3:prm.nx+3] ) # XXX check bounds

    return ajx



