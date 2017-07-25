import numpy as np

def efield(ex, ajx, prm):

    #amperes law E_n+1 = E_n - J
    ex[3:prm.nx+3] = ex[3:prm.nx+3] - ajx[3:prm.nx+3]

    #wrap
    ex[0:2]               = ex[prm.nx:prm.nx+2]
    ex[prm.nx+3:prm.nx+4] = ex[3:4]

    #shift
    fex = np.zeros(prm.nx + 6)
    for ii in range(3, prm.nx+3):
        fex[ii] = (ex[ii] + ex[ii-1])/2.0

    #wrap
    fex[0:2]               = fex[prm.nx:prm.nx+2]
    fex[prm.nx+3:prm.nx+4] = fex[3:4]

    return ex, fex



def efield_f(fex, ajx, prm):
    fex[3:prm.nx+3] = fex[3:prm.nx] - ajx[3:prm.nx+3]

    #wrap
    fex[0:2]               = fex[prm.nx:prm.nx+2]
    fex[prm.nx+3:prm.nx+4] = fex[3:4]

    return fex


