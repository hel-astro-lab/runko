import numpy as np
import h5py
import glob
import re
import os

from configSetup import Configuration

from combine_files import natural_keys
from combine_files import atoi
from combine_files import combine_tiles


# merge field file
def merge_field_nodes(fdir, lap, conf):

    fstr = fdir + '/fields-*_'+str(lap)+'.h5'
    print(fstr)
    files = glob.glob(fstr) 
    files.sort(key=natural_keys)

    print(files)

    fname_merged = fdir + '/fields_'+str(lap)+'.h5'
    f5 = h5py.File(fname_merged,'w')

    fvars = [
            'bx','by','bz',
            'ex','ey','ez',
            'jx','jy','jz',
            'rho'
            ]

    for fvar in fvars:
        print("processing {}".format(fvar))

        arr = combine_tiles(files[0], fvar, conf, isp=None)
        for ff in files[1:]:
            print("  file: {}".format(ff))
            arrN = combine_tiles(ff, fvar, conf, isp=None)
            arr += arrN

        dset = f5.create_dataset(fvar, data=arr)

    f5.close()


# merge field file
def merge_analysis_nodes(fdir, lap, isps, conf):

    fstr = fdir + '/analysis-*_'+str(lap)+'.h5'
    print(fstr)
    files = glob.glob(fstr) 
    files.sort(key=natural_keys)

    print(files)

    fname_merged = fdir + '/analysis_'+str(lap)+'.h5'
    f5 = h5py.File(fname_merged,'w')

    fvars = [
            'rho',
            'edens',
            #'temp',
            #'Vx',
            #'Vy',
            #'Vz',
            #'momx',
            #'momy',
            #'momz',
            #'pressx',
            #'pressy',
            #'pressz',
            #'shearxy',
            #'shearxz',
            #'shearyz',
            ]

    for fvar in fvars:
        print("processing {}".format(fvar))

        for isp in isps:
            print(" isp {}".format(isp))
            arr = combine_tiles(files[0], fvar, conf, isp=isp)
            for ff in files[1:]:
                print("  file: {}".format(ff))
                arrN = combine_tiles(ff, fvar, conf, isp=None)
                arr += arrN

            dset = f5.create_dataset(fvar+'_'+str(isp), data=arr)

    f5.close()



if __name__ == "__main__":


    #conf = Configuration('turb100hot.ini')
    fdir = "out100hot"

    lap = 100
    conf = Configuration('config-test.ini')
    fdir = 'out'

    merge_field_nodes(fdir, lap, conf)
    merge_analysis_nodes(fdir, lap, [0,1], conf)




