import numpy as np
import os, sys
import copy
import matplotlib.pyplot as plt
import palettable as pal
from matplotlib import cm

from visualize import *
import vmesh


#Profiling / Timing
sys.path.insert(0, '../tools')
from timer import Timer


#physical "real" distribution to compare against
def physical_vel(x,y,z):

    mux = 1.0
    muy = 2.0
    muz = 3.0
    sigmax = 2.0
    sigmay = 3.0
    sigmaz = 4.0

    vx = np.exp(-(x-mux)**2 / sigmax**2 )
    vy = np.exp(-(y-muy)**2 / sigmay**2 )
    #vz = np.exp(-(z-muz)**2 / sigmaz**2 )
    vz = 1.0

    return vx*vy*vz
    #return 0.5


def populate_mesh( mesh ):

    for k in range(mesh.Nblocks[2]):
        for j in range(mesh.Nblocks[1]):
            for i in range(mesh.Nblocks[0]):
                cid = mesh.get_block_ID([i,j,k])
                (x,y,z) = mesh.get_center( cid )

                fval = physical_vel(x,y,z)
                mesh[i,j,k] = [fval, fval, fval, fval]
                #print "({},{},{}) = {}".format(i,j,k,fval)


class Params:
    mins = None
    maxs = None
    lens = None




if __name__ == "__main__":


    ################################################## 
    # set up plotting and figure
    plt.fig = plt.figure(1, figsize=(8,4))
    plt.rc('font', family='serif', size=12)
    plt.rc('xtick')
    plt.rc('ytick')
    
    gs = plt.GridSpec(1, 2)
    gs.update(hspace = 0.5)
    
    axs = []
    axs.append( plt.subplot(gs[0]) )
    axs.append( plt.subplot(gs[1]) )

    
    #cid = 1
    #for zi in range(10):
    #    for yi in range(10):
    #        for xi in range(10):
    #            GID = zi*(10*10) + yi*10 + xi + 1
    #            print "({},{},{}) cid: {} GID: {}".format(xi,yi,zi,cid,GID)
    #            cid += 1
    

    ################################################## 
    # set-up grid
    # xy
    params = Params()
    params.mins = [ -20.0, -20.0, -20.0 ]
    params.maxs = [  20.0,  20.0,  20.0 ]

    mesh = vmesh.vMesh()
    mesh.Nblocks = [50, 30, 5]
    #mesh.Nblocks = [100, 6, 2]
    mesh.zFill(params.mins, params.maxs)

    populate_mesh( mesh )

    print "num of blocks:", mesh.number_of_blocks
    #cid = mesh.get_block_ID( [5,5,5] )
    #print "cid:",cid
    #print "cen:",mesh.get_center(cid)
    #print "val:",mesh[5,5,5]
    print "Memory usage: {} Mb ({} cells)".format(mesh.sizeInBytes()/1e6, mesh.number_of_blocks)


    #mesh.clip()
    print "num of blocks:", mesh.number_of_blocks

    #visualize_mesh(axs[0], mesh, params)
    #visualize_data(axs[1], mesh, params)


    print "Memory usage: {} Mb ({} cells)".format(mesh.sizeInBytes()/1e6, mesh.number_of_blocks)


    ##################################################
    # test bundle math
    print "testing bundle..."
    vbundle = mesh.get_bundle(0, 0, 0)
    print vbundle.getGrid()
    print vbundle.getPencil()

    #mesh.add_bundle(0, 7, 0, vbundle)
    #vbundle = mesh.get_bundle(0, 20, 0)

    #visualize_data2(axs[0], mesh, params)
    #plt.savefig("vlasov_0001.png")


    ##################################################
    #test propagation
    vbundle = mesh.get_bundle(0, 20, 0)
    #intp = vmesh.BundleInterpolator2nd()
    #intp.setBundle(vbundle)



    vbundle = mesh.get_bundle(2, 20, 0)
    mesh.add_bundle(2, 10, 0, vbundle)


    visualize_data2(axs[0], mesh, params)
    visualize_data3(axs[1], mesh, params)
    plt.savefig("vlasov0.png")
    sys.exit()
    ##################################################




    vsol = vmesh.vSolver()
    vsol.setMesh(mesh)
    #intp = vmesh.BundleInterpolator2nd()
    intp = vmesh.BundleInterpolator4th()
    vsol.setInterpolator(intp)


    timer = Timer(["vsol"])
    timer.start("vsol")

    for lap in range(100000):
        vsol.solve()
        #vsol.vmesh.clip()
        timer.lap("vsol")

        if (lap % 200 == 0): 
            print("lap : {}").format(lap)
            visualize_data2(axs[0], vsol.vmesh, params)
            stri = str(lap).rjust(4, '0')

            #measure temperature
            #bundle = mesh.get_bundle(0, 25, 2)

            #axs[1].plot( bundle.getGrid(), bundle.getPencil() )

            visualize_data(axs[1], vsol.vmesh, params)
            plt.savefig("out/vlasov_"+stri+".png")

    timer.stats("vsol")






