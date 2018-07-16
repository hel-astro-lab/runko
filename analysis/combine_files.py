import numpy as np
import h5py
import glob
import re

from configSetup import Configuration


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]



# list & organize output files
def get_file_list(fdir, fname):
    fstr = fdir + fname + "*.h5"
    files = glob.glob(fstr) 

    files.sort(key=natural_keys)

    return files



# combine different time snapshots together
def combine_files(fdir, 
                  fname, 
                  var, 
                  conf, 
                  isp=None, 
                  func=None):

    fstr = fdir + fname + "*.h5"
    files = glob.glob(fstr) 

    files.sort(key=natural_keys)
    #print(files)

    ret_arr = []
    
    for ff in files:
        arr = combine_tiles(ff, var, conf, isp=isp)
        
        if not(func == None):
            arr = func(arr)
        ret_arr.append( arr )

    return np.array(ret_arr)



# combine tiles inside node together into one array
def combine_tiles(ff, fvar, conf, isp=None ):

    arr = np.zeros((conf.Nx*conf.NxMesh, 
                    conf.Ny*conf.NyMesh, 
                    conf.Nz*conf.NzMesh))

    f = h5py.File(ff,'r')
    for dset in f:

        if not(isp==None):
            if not(f[dset]['ispcs'].value == isp):
                continue

        i = f[dset]['i'].value
        j = f[dset]['j'].value
        k = f[dset]['k'].value

        NxMesh = f[dset]['Nx'].value
        NyMesh = f[dset]['Ny'].value
        NzMesh = f[dset]['Nz'].value
    
        ii = int( i*NxMesh )
        jj = int( j*NyMesh )
        kk = int( k*NzMesh )

        tile = f[dset][fvar][()]
        tile = np.reshape(tile, (NzMesh, NyMesh, NxMesh))

        for s in range(NzMesh):
            for r in range(NyMesh):
                for q in range(NxMesh):
                    arr[ii+q, jj+r, kk+s] = tile[s,r,q]




    return arr


##################################################

if __name__ == "__main__":

    conf = Configuration('../pic/config-weibel.ini') 
    fdir = "../pic/weibel/out/"
    #fname = "field"
    #var = "ex"
    
    fname = "analysis"
    var = "mgamma"
    
    #arrs = combine_files(fdir, fname, var, conf)
    #arrs = combine_files(fdir, fname, var, conf, func=np.max)
    arrs = combine_files(fdir, fname, var, conf, isp=0)
    
    print(arrs)
    print(arrs.shape)
    print(len(arrs))
    
    #ftest = "../pic/weibel/out/fields-0_0.h5"
    #arr = combine_tiles(ftest, var, conf)
    #print(arr[:,:,0])









