import numpy as np

import sys, os

import corgi
import pyplasma as plasma

np.random.seed(0)



def spatialLoc(node, Ncoords, Mcoords, conf):

    #node coordinates
    i, j    = Ncoords 
    Nx      = conf.Nx
    Ny      = conf.Ny

    #mesh coordinates
    l, m, n = Mcoords 
    NxMesh = conf.NxMesh
    NyMesh = conf.NyMesh
    NzMesh = conf.NzMesh

    #grid spacing
    xmin = node.getXmin()
    ymin = node.getYmin()

    dx = conf.dx
    dy = conf.dy
    dz = conf.dz


    #calculate coordinate extent
    x = xmin + i*NxMesh*dx + l*dx
    y = ymin + j*NyMesh*dy + m*dy
    z = 0.0                + n*dz

    return (x, y, z)



def createEmptyVelocityMesh(conf):
    vmesh = plasma.AdaptiveMesh3D()

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


def fillMesh(vmesh, ffunc, xloc, ispcs, conf):

    # standard flat fill on 0th rfl level
    nx, ny, nz = vmesh.get_size(0)
    for r in range(nx):
        for s in range(ny):
            for t in range(nz):
                uloc = vmesh.get_center([r,s,t], 0)

                val = ffunc(xloc, uloc, ispcs, conf)
                vmesh[r,s,t, 0] =  val #ref lvl 0


    if conf.refinement_level < 1:
        if conf.clip:
            vmesh.clip_cells(conf.clipThreshold)
        return 


    ###################################################
    # adaptivity

    adapter = plasma.Adapter();

    sweep = 1
    while(True):
        #print("-------round {}-----------".format(sweep))
        adapter.check(vmesh)
        adapter.refine(vmesh)

        #print("cells to refine: {}".format( len(adapter.cells_to_refine)))
        for cid in adapter.cells_created:
            rfl = vmesh.get_refinement_level(cid)
            indx = vmesh.get_indices(cid)
            uloc = vmesh.get_center(indx, rfl)

            val = ffunc(xloc, uloc, ispcs, conf)
            vmesh[indx[0], indx[1], indx[2], rfl] = val

        adapter.unrefine(vmesh)
        #print("cells to be removed: {}".format( len(adapter.cells_removed)))

        sweep += 1
        if sweep > conf.refinement_level: break

    if conf.clip:
        vmesh.clip_cells(conf.clipThreshold)

    return 



#inject plasma into cells
def inject(node, ffunc, conf):

    #loop over all *local* cells
    for i in range(node.getNx()):
        for j in range(node.getNy()):
            #if n.getMpiGrid(i,j) == n.rank:
            if True:
                print("creating ({},{})".format(i,j))

                #get cell & its content
                cid    = node.cellId(i,j)
                c      = node.getCellPtr(cid) #get cell ptr

                # loop over species
                species = []
                for ispcs in range(2):
                    block = plasma.PlasmaBlock(conf.NxMesh, conf.NyMesh, conf.NzMesh)
                    
                    #set q/m
                    if ispcs == 0:
                        block.qm = conf.me
                    elif ispcs == 1:
                        block.qm = conf.mi


                    for n in range(conf.NzMesh):
                        for m in range(conf.NyMesh):
                            for l in range(conf.NxMesh):
                                #print(" sub mesh: ({},{},{})".format(l,m,n))

                                xloc = spatialLoc(node, (i,j), (l,m,n), conf)

                                vmesh = createEmptyVelocityMesh(conf)
                                fillMesh(vmesh,
                                         ffunc,
                                         xloc,
                                         ispcs,
                                         conf)
                                #vmesh.maximum_refinement_level = conf.refinement_level


                                block[l,m,n] = vmesh
                    species.append(block)

                c.insertInitialSpecies(species)

