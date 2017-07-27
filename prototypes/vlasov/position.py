import numpy as np


def position_linear(ff, vx, prm):
    ajxs = np.zeros( (prm.nx+6, prm.ns) )
    flux = np.zeros( (prm.nv+6, prm.nx+6) )

    for kk in range(prm.ns):
        aa = vx[:, kk] / prm.dx * prm.dt
        fa = np.floor(aa)

        for ii in range(2, prm.nx+3):
            iss = np.ones(prm.nv+6)*ii - fa
            flux[:, ii] = aa[:] * np.diag(ff[:, iss.astype(int), kk])

        ff[:, 3:prm.nx+3, kk] = ff[:, 3:prm.nx+3, kk] \
                - (flux[:, 3:prm.nx+3] - flux[:, 2:prm.nx+2])
                
    #wrap
    ff[:, 0:2, :]               = ff[:, prm.nx:prm.nx+2, :]
    ff[:, prm.nx+3:prm.nx+6, :] = ff[:, 3:6, :]
                    
    return ff



