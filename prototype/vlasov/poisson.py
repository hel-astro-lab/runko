import numpy as np


def poisson(ex, rho, prm):

    #XXX +3 or +2
    for ii in range(3, prm.nx+3):
        ex[ii] = ex[ii-1] + rho[ii]

    #wrap boundaries
    ex[0:2]               = ex[prm.nx:prm.nx+2]
    ex[prm.nx+3:prm.nx+4] = ex[3:4]

    ex0 = np.sum( ex[3:prm.nx+3] ) / prm.nx
    ex -= ex0

    #relocate
    fex = np.zeros(prm.nx+6) #already defined
    for ii in range(3, prm.nx+2):
        fex[ii] = ( ex[ii] + ex[ii-1] )/2.0

    fex[0:2]               = fex[prm.nx:prm.nx+2]
    fex[prm.nx+3:prm.nx+6] = fex[3:6]

    return ex, fex
