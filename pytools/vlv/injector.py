# -*- coding: utf-8 -*-

import numpy as np
import sys

import pyrunko


def ind2loc(grid, gridI, tileI, conf):

    # grid coordinates
    i, j = gridI
    Nx = conf.Nx
    Ny = conf.Ny

    # tile coordinates
    l, m, n = tileI
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    # grid spacing; start point + step
    xmin = conf.xmin
    ymin = conf.ymin
    zmin = 0.0

    dx = 1.0  # conf.dx
    dy = 1.0  # conf.dy
    dz = 1.0  # conf.dz

    # calculate coordinate extent
    x = xmin + i * (NxMesh) * dx + l * dx
    y = ymin + j * (NyMesh) * dy + m * dy
    z = zmin + 0 * (NzMesh) * dz + n * dz

    return [x, y, z]



def empty_filler(xloc, uloc, ispcs, conf):
    return 0.0

def createEmptyVelocityMesh(conf):
    vmesh = pyrunko.tools.AdaptiveMesh3D()

    dx = (conf.vxmax - conf.vxmin)/(conf.Nvx)
    dy = (conf.vymax - conf.vymin)/(conf.Nvy)
    dz = (conf.vzmax - conf.vzmin)/(conf.Nvz)

    if (conf.Nvy == 1) and (conf.Nvz == 1): #1D switch
        vmesh.resize( [conf.Nvx,  2,  2 ])
        vmesh.set_min([conf.vxmin,    -0.5, -0.5])
        vmesh.set_max([conf.vxmax-dx,  0.5,  0.5])

        vmesh.top_refinement_level = conf.refinement_level
        return vmesh

    else:
        vmesh.resize( [conf.Nvx,      conf.Nvy,      conf.Nvz ])
        vmesh.set_min([conf.vxmin,    conf.vymin,    conf.vzmin])
        vmesh.set_max([conf.vxmax-dx, conf.vymax-dy, conf.vzmax-dz])

        vmesh.top_refinement_level = conf.refinement_level
        return vmesh


def fillMesh(
        vmesh, 
        ffunc, 
        xloc, 
        ispcs, 
        conf, 
        preclip = lambda a,b,c,d : False
        ):

    nx, ny, nz = vmesh.get_size(0)

    #collapse fill to 1V 
    if (ny < 3) and (nz < 3):
        ny = 1
        nz = 1

    # standard flat fill on 0th rfl level
    for t in range(nz):
        for s in range(ny):
            for r in range(nx):
                uloc = vmesh.get_center([r,s,t], 0)

                if preclip(xloc, uloc, ispcs, conf): continue

                val = ffunc(xloc, uloc, ispcs, conf)
                vmesh[r,s,t, 0] =  val #ref lvl 0


    if conf.refinement_level < 1:
        if conf.clip:
            vmesh.clip_cells(conf.clipThreshold)
        return 


    ###################################################
    # adaptivity

    adapter = pyrunko.Adapter();

    sweep = 1
    while(True):
        #print("-------round {}-----------".format(sweep))
        adapter.check(vmesh)
        adapter.refine(vmesh)

        #print("tiles to refine: {}".format( len(adapter.tiles_to_refine)))
        for cid in adapter.tiles_created:
            rfl = vmesh.get_refinement_level(cid)
            indx = vmesh.get_indices(cid)
            uloc = vmesh.get_center(indx, rfl)

            val = ffunc(xloc, uloc, ispcs, conf)
            vmesh[indx[0], indx[1], indx[2], rfl] = val

        adapter.unrefine(vmesh)
        #print("tiles to be removed: {}".format( len(adapter.tiles_removed)))

        sweep += 1
        if sweep > conf.refinement_level: break

    if conf.clip:
        vmesh.clip_tiles(conf.clipThreshold)

    return 


def inject_internal(i,j,
                    grid, 
                    ffunc, 
                    conf,
                    preclip = lambda a,b,c,d : False
                    ):

    print("creating parallel ({},{})".format(i,j))
    
    #get tile & its content
    cid    = grid.id(i)
    c      = grid.get_tile(cid) #get tile ptr
    
    # loop over species
    species = []
    for ispcs in range(conf.Nspecies):
        block = pyrunko.vlv.PlasmaBlock(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        
        #set q/m
        if ispcs == 0:
            block.qm = conf.me
        elif ispcs == 1:
            block.qm = conf.mi
        elif ispcs == 2:
            block.qm = conf.me
        elif ispcs == 3:
            block.qm = conf.mi
    
    
        for n in range(conf.NzMesh):
            for m in range(conf.NyMesh):
                for l in range(conf.NxMesh):
                    #print(" sub mesh: ({},{},{})".format(l,m,n))
    
                    xloc = ind2loc(grid, (i,j), (l,m,n), conf)
    
                    vmesh = createEmptyVelocityMesh(conf)
                    fillMesh(vmesh,
                             ffunc,
                             xloc,
                             ispcs,
                             conf,
                             preclip=preclip,
                             )
                    #vmesh.maximum_refinement_level = conf.refinement_level
    
                    block[l,m,n] = vmesh
        species.append(block)
    
    c.insert_initial_species(species)


#inject plasma into tiles
def inject_parallel(
        grid, 
        ffunc, 
        conf,
        preclip = lambda a,b,c,d : False
        ):

    #multiprocessing 
    pool = Pool() 
    nxnynz = [(i,j) for i in range(grid.get_Nx()) for j in range(grid.get_Ny()) ]
    print("pool for injector:", nxnynz)
    #pool.map(inject_internal, nxnynz)
    pool.map(partial(inject_internal,
                        grid=node,
                        ffunc=ffunc,
                        conf=conf,
                        preclip=preclip
                        ),
                    nxnynz)


#inject plasma into tiles
def inject(
        grid, 
        ffunc, 
        conf,
        preclip = lambda a,b,c,d : False,
        empty = False,
        ):

    # setup toolbar
    toolbar_width = conf.Nx*conf.Ny
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

    #loop over all *local* tiles
    for i in range(grid.get_Nx()):
        for j in range(grid.get_Ny()):

            #if n.get_mpi_grid(i,j) == n.rank:
            if True:
                #print("creating ({},{})".format(i,j))
                sys.stdout.write("-")
                sys.stdout.flush()

                #get tile & its content
                cid    = grid.id(i)
                c      = grid.get_tile(cid) #get tile ptr

                # loop over species
                species = []
                for ispcs in range(conf.Nspecies):
                    block = pyrunko.vlv.PlasmaBlock(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                    
                    #set q/m
                    if ispcs == 0:
                        block.qm = conf.me
                    elif ispcs == 1:
                        block.qm = conf.mi
                    elif ispcs == 2:
                        block.qm = conf.me
                    elif ispcs == 3:
                        block.qm = conf.mi


                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            for l in range(conf.NxMesh):
                                #print(" sub mesh: ({},{},{})".format(l,m,n))

                                xloc = ind2loc(grid, (i,j), (l,m,n), conf)

                                vmesh = createEmptyVelocityMesh(conf)
                                if not(empty):
                                    fillMesh(vmesh,
                                             ffunc,
                                             xloc,
                                             ispcs,
                                             conf,
                                             preclip=preclip,
                                             )
                                #vmesh.maximum_refinement_level = conf.refinement_level

                                block[l,m,n] = vmesh
                    species.append(block)

                c.insert_initial_species(species)

    sys.stdout.write("\n") #finish toolbar




