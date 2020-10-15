# -*- coding: utf-8 -*- 


# dig out indices from inside the tile
def get_index(tile, conf):
    ind = tile.index
    if conf.threeD:
        i,j,k = ind
    elif conf.twoD:
        i,j = ind
        k = 0

    return i,j,k


class Stagger:

    staggers = {
            'ex':[1,0,0],
            'ey':[0,1,0],
            'ez':[0,0,1],

            'bx':[0,1,1],
            'by':[1,0,1],
            'bz':[1,1,0],

            'jx':[1,0,0],
            'jy':[0,1,0],
            'jz':[0,0,1],

            'rh':[0,0,0],

               }

    # compute transformation indices for going from
    # variable at loc1 staggering to loc2 staggering
    def x2y(self, loc1, loc2):
        if loc2 == 'no':
            return [0.0, 0.0, 0.0]

        offs1 = self.staggers[loc1]
        offs2 = self.staggers[loc2]  
        
        ret = [
                -0.5*(offs2[0] - offs1[0]),
                -0.5*(offs2[1] - offs1[1]),
                -0.5*(offs2[2] - offs1[2]),
              ]

        return ret

    def at(self, stg, stg0='rh'):
        #if stg == 'no' or stg0 == 'no': #no staggering at all; just dummy copy of self 
        #    ret = Stagger(self.x, self.y ,self.z)

        offs = self.x2y(stg0, stg)

        ret = Stagger(self.x, self.y ,self.z)
        ret.x += offs[0]
        ret.y += offs[1]
        ret.z += offs[2]
        return ret


    def __init__(self, x,y,z):
        self.x = x
        self.y = y
        self.z = z

    #overloaded operators

    #+ float/int
    #- float/int
    #* float/int
    #/ float/int
    #** float/int
    #+ stg
    #- stg

    def __add__(self, o):
        # always upcast return object to Stagger class
        ret = Stagger(self.x, self.y ,self.z)

        #add component wise
        if isinstance(o, self.__class__):
            ret.x += o.x
            ret.y += o.y
            ret.z += o.z

            return ret
        else:
            raise TypeError("unsupported operand type(s) for +: '{}' and '{}'").format(self.__class__, type(o)) 

    def __sub__(self, o):
        # always upcast return object to Stagger class
        ret = Stagger(self.x, self.y ,self.z)

        #subtract component wise
        if isinstance(o, self.__class__):
            ret.x -= o.x
            ret.y -= o.y
            ret.z -= o.z

            return ret
        else:
            raise TypeError("unsupported operand type(s) for -: '{}' and '{}'").format(self.__class__, type(o)) 




def ind2loc(gridI, tileI, conf):

    # grid coordinates
    i, j, k = gridI
    Nx = conf.Nx
    Ny = conf.Ny
    Nz = conf.Nz

    # tile coordinates
    l, m, n = tileI
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    # grid spacing; start point + step
    xmin = conf.xmin
    ymin = conf.ymin
    zmin = conf.zmin

    dx = 1.0  # conf.dx
    dy = 1.0  # conf.dy
    dz = 1.0  # conf.dz

    # calculate coordinate extent
    x = xmin + i * (NxMesh) * dx + l * dx
    y = ymin + j * (NyMesh) * dy + m * dy
    z = zmin + k * (NzMesh) * dz + n * dz

    return [x, y, z]





















