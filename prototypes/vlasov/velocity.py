import numpy as np


def velocity_linear(ff, ex, prm):
    flux = np.zeros( (prm.nv+6, prm.nx+6) )

    jj = np.arange(prm.nv+6)
    for kk in range(prm.ns):
        aa = ex[:] * prm.qm[kk]/prm.dv[kk]*prm.dt
        fa = np.floor(aa).astype(int)

        for ii in range(prm.nx+5):
            js = jj - fa[ii] - 1
            flux[jj, ii] = aa[ii] * ff[js, ii, kk]

        ff[1:prm.nv+5, :, kk] = ff[1:prm.nv+5, :, kk] \
                - (flux[1:prm.nv+5, :] - flux[0:prm.nv+4, :])

    return ff
