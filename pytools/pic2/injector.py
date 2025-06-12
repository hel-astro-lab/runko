# -*- coding: utf-8 -*-

# dummy function that returns always 1 (for weight initialization)
def unit_w(xloc, ispcs, conf):
    return 1.0


# inject plasma into (individual) cells
def inject(grid,
           vel_func,
           den_func,
           conf,
           align_species=True,
           w_func=unit_w,
           ):
    print("warning: pytools.pic2.inject called but not implemented")
    return 42
