import numpy as np
from initial import initial
from charge import charge
from poisson import poisson

from position import position_linear
from velocity import velocity_linear
from current import current
from efield import efield



#initialize
#-------------------------------------------------- 
#load configuration
import twostream as prm

#hdiag = diagnostics_init(prm)



#initial(prm)
#ff, gx, gv, ex, fex, ajx, xx, vx, kx, kv, ifdiag = initial(prm)
ff, gx, gv, ex, fex, ajx, xx, vx, kx, kv = initial(prm)

#initial step
rho = charge(ff, prm)
ex,fex = poisson(ex, rho, prm)
ff = position_linear(ff, vx, prm)
ajx = current(ff, vx, prm)


#conservative schemes 
# 0 linear
# 1 2nd order
# 2 4th order
# 6 CPIC4
ex, fex = efield(ex, ajx, prm)

#non-conservative
# 3 - cubic spline
# 5 - CIP3 
#rho = charge(ff, prm)
#ex, fex = poisson(ex, rho, prm)

#3rd possible case ???
#fex = efield_f(fex, ajx, prm)


#-------------------------------------------------- 
# main loop

jtime = 0
for jtime in range(prm.ntime):
    print "-----------", jtime, "----------"

    ff  = velocity_linear(ff, fex, prm)
    ff  = position_linear(ff, vx, prm)

    ajx = current(ff, vx, prm)
    rho = charge(ff, prm)

    ex, fex = efield(ex, ajx, prm)
    











