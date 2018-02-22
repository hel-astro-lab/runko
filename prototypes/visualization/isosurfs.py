from __future__ import print_function

import numpy as np
from tvtk.api import tvtk
from mayavi.scripts import mayavi2

#import pyplasma as plasma
#import pyplasmaDev as pdev

from injector import createEmptyVelocityMesh, fillMesh
from scipy.stats import multivariate_normal


class Conf:

    outdir = "out"

    #-------------------------------------------------- 
    # velocity space parameters
    vxmin = -10.0
    vymin = -10.0
    vzmin = -10.0

    vxmax =  10.0
    vymax =  10.0
    vzmax =  10.0

    Nvx = 25
    Nvy = 25
    Nvz = 25

    #vmesh refinement
    refinement_level = 0
    clip = True
    clipThreshold = 1.0e-5



def filler(xloc, uloc, ispcs, conf):

    delgam = np.sqrt(1.0)

    mean = [ 0.0, 0.0, 0.0]
    cov  = np.zeros((3,3))
    cov[0,0] = 1.0 
    cov[1,1] = 1.0 
    cov[2,2] = 6.0

    f = multivariate_normal.pdf(uloc, mean, cov)

    mean = [ 1.0, 1.0, 1.0]
    cov  = np.zeros((3,3))
    cov[0,0] = 4.0 
    cov[1,1] = 3.0 
    cov[2,2] = 2.0
    f += multivariate_normal.pdf(uloc, mean, cov)

    return f


# Strip AMR velocity mesh into something that VTK understands
def stripAmrMesh(vmesh):

    cids = vmesh.get_cells(True)
    Ncells = len(cids)

    Ncells = 10

    points = np.zeros((Ncells, 3))
    conns = np.zeros((Ncells, 6))
    vals   = np.zeros(Ncells)

    for q,cid in enumerate(cids):
        if q > 9:
            break
        print(q)

        i,j,k = vmesh.get_indices(cid)
        rfl   = vmesh.get_refinement_level(cid)
        xloc = vmesh.get_center([i,j,k], rfl)
        dxs  = vmesh.get_length(rfl)

        #ip = 0
        #for ix in [-1.0, 1.0]:
        #    for iy in [-1.0, 1.0]:
        #        for iz in [-1.0, 1.0]:
        #            xx = xloc[0] + ix*dxs[0]/2
        #            yy = xloc[1] + iy*dxs[1]/2
        #            zz = xloc[2] + iz*dxs[2]/2

        #            points[q+ip, :] = [xx, yy, zz]
        #            #conns[q, 
        #            ip += 1

        points[q, :] = xloc
        vals[q] = vmesh[i,j,k,rfl]

    return points, conns, vals


def denseMatrix(vmesh):
    cids = vmesh.get_cells(True)
    Ncells = len(cids)
    Nx, Ny, Nz = vmesh.get_size(0)
    data = np.zeros((Nx, Ny, Nz))

    for q,cid in enumerate(cids):
        i,j,k = vmesh.get_indices(cid)
        rfl   = vmesh.get_refinement_level(cid)
        xloc  = vmesh.get_center([i,j,k], rfl)

        #dxs  = vmesh.get_length(rfl)
        #ip = 0
        #for ix in [-1.0, 1.0]:
        #    for iy in [-1.0, 1.0]:
        #        for iz in [-1.0, 1.0]:
        #            xx = xloc[0] + ix*dxs[0]/2
        #            yy = xloc[1] + iy*dxs[1]/2
        #            zz = xloc[2] + iz*dxs[2]/2

        #            points[q+ip, :] = [xx, yy, zz]
        #            #conns[q, 
        #            ip += 1

        val = vmesh[i,j,k,rfl]

        data[i,j,k] = val

    return data

def image_data(vmesh):

    data = denseMatrix(vmesh)
    i = tvtk.ImageData(spacing=(1, 1, 1), origin=(0, 0, 0))
    i.point_data.scalars = data.ravel()
    i.point_data.scalars.name = 'scalars'
    i.dimensions = data.shape

    return i



def createVtkObject(vmesh):
    print("Creating TVTK object from AMR mesh")

    #points = array([[0,1.2,0.6], [1,0,0], [0,1,0], [1,1,1], # tetra
    #                [1,0,-0.5], [2,0,0], [2,1.5,0], [0,1,0],
    #                [1,0,0], [1.5,-0.2,1], [1.6,1,1.5], [1,1,1], # Hex
    #                ], 'f')

    ## The cells
    #cells = array([4, 0, 1, 2, 3, # tetra
    #               8, 4, 5, 6, 7, 8, 9, 10, 11 # hex
    #               ])

    ## The offsets for the cells, i.e. the indices where the cells start.
    #offset = array([0, 5])
    #tetra_type = tvtk.Tetra().cell_type # VTK_TETRA == 10
    #hex_type = tvtk.Hexahedron().cell_type # VTK_HEXAHEDRON == 12
    #cell_types = array([tetra_type, hex_type])

    ## Create the array of cells unambiguously.
    #cell_array = tvtk.CellArray()
    #cell_array.set_cells(2, cells)

    ## Now create the UG
    #ug = tvtk.UnstructuredGrid(points=points)
    #ug.set_cells(cell_types, offset, cell_array)
    #scalars = random.random(points.shape[0])
    #ug.point_data.scalars = scalars
    #ug.point_data.scalars.name = 'scalars'


    #analysator
    #points = nodesAndKeys[0]
    #tets = nodesAndKeys[1]

    # [x,y,z] list of nodes
    #points = np.array([[0,0,0],
    #                   [1,0,0],
    #                   [0,1,0],
    #                   [1,1,0],
    #                   [0,0,1],
    #                   [1,0,1],
    #                   [1,1,1],
    #                   ], 'f')

    #cells = np.array([[0,1,2,3],
    #                  [4,5,6,7],
    #                  [4,5,6,7],


    points, conns, avgs = stripAmrMesh(vmesh)

    print(points)
    print(conns)

    #set cells & grid
    vox_type = tvtk.Voxel().cell_type#VTK_VOXEL
    ug=tvtk.UnstructuredGrid(points=points)
    ug.set_cells(vox_type, conns)


    #points = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1], # tets
    #                [1,0,0], [2,0,0], [1,1,0], [1,0,1],
    #                [2,0,0], [3,0,0], [2,1,0], [2,0,1],
    #                ], 'f')
    #tets = np.array([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]])
    #tet_type = tvtk.Tetra().cell_type
    #ug = tvtk.UnstructuredGrid(points=points)
    #ug.set_cells(tet_type, tets)

    #set scalar values
    #avgs = np.array([10.0])
    values=np.ravel(avgs)
    ug.cell_data.scalars=values
    ug.cell_data.scalars.name='dens'


    #d = mayavi.mlab.pipeline.add_dataset(ug)

    #if iso_surface == False:
    #    iso = mayavi.mlab.pipeline.surface(d)
    #else:
    #ptdata = mayavi.mlab.pipeline.cell_to_point_data(d)
    #iso = mayavi.mlab.pipeline.iso_surface(ptdata, contours=[1e-15,1e-14,1e-12], opacity=0.3)


    return ug



# create AMR mesh
conf = Conf()
vmesh = createEmptyVelocityMesh(conf)
fillMesh(vmesh, filler, [0.0, 0.0, 0.0], 0, conf)


# create VTK objects
#ug = createVtkObject(vmesh)
ug = image_data(vmesh)



def view_mlab():
    engine = mayavi.mlab.get_engine()

    from mayavi.modules.axes import Axes 
    axes = Axes()
    axes.name = 'Axes'
    axes.axes.fly_mode = 'none'
    axes.axes.number_of_labels = 8
    axes.axes.font_factor = 0.5
    #module_manager = self.__module_manager()
    # Add the label / marker:
    engine.add_filter( axes )
    from mayavi.modules.outline import Outline
    outline = Outline()
    outline.name = 'Outline'
    engine.add_filter( outline )

    mayavi.mlab.show()




@mayavi2.standalone
def view():

    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi.modules.outline import Outline
    from mayavi.modules.surface import Surface
    from mayavi.modules.vectors import Vectors
    from mayavi.modules.api import IsoSurface

    scene = mayavi.new_scene()
    #scene.background = "black"


    
    # The single type one
    src = VTKDataSource(data = ug)
    mayavi.add_source(src)
    mayavi.add_module(Outline())
    #mayavi.add_module(Surface())
    #mayavi.add_module(Vectors())

    #translucent isosurfaces
    iso = IsoSurface()
    mayavi.add_module(iso)
    iso.module_manager.scalar_lut_manager.lut_mode = "hot"

    iso.contour.contours = np.linspace(0.0, 0.03, 30).tolist()
    iso.actor.property.opacity = 0.3


    #iso = mayavi.mlab.pipeline.iso_surface(ug, contours=[1e-15,1e-14,1e-12], opacity=0.3)
    #from mayavi import mlab
    #mlab.contour3d(ug, opacity=0.3)



if __name__ == '__main__':
    view()


