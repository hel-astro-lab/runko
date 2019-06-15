import numpy as np
import h5py
import glob
import re
import os

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
    fstr = fdir + fname + '_*.h5'
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

    fstr = fdir + fname + '-*_*.h5'
    files = glob.glob(fstr) 
    files.sort(key=natural_keys)

    #last file
    #print(files)
    #print(files[-1])

    if fname == 'analysis':
        fname_prepro = fdir + "processed-" + fname + "-" + var + "_" + str(isp) + ".h5"
    elif fname == 'fields':
        fname_prepro = fdir + "processed-" + fname + "-" + var + ".h5"

    #print(fname_prepro)

    #print( re.split('(\d+)', files[-1].strip()) )
    last_file_num = atoi( re.split('(\d+)', files[-1].strip())[-4] )
    print("...newest file is: ", last_file_num)

    if os.path.isfile(fname_prepro):
        print("...reading from preprocesses file: ", fname_prepro)
        f5 = h5py.File(fname_prepro,'r')
        dset = f5['processed']
        newest_snapshot = dset.attrs['current_head']
        print("...latest processed file is ", newest_snapshot)

        #use what is in the file and read the rest, then update the tail 
        if newest_snapshot < last_file_num:
            #print(" adding files from {} to {}".format(newest_snapshot, last_file_num))

            # skim file list so that only newest additions are left
            # this way we do not reprocess stuff

            i = 0
            for ff in files:
                file_num = atoi( re.split('(\d+)', ff.strip())[-4] )
                if file_num > newest_snapshot:
                    break
                else:
                    i += 1
            #print(" using only files")
            #print(files[i:])
            files = files[i:]

            #TODO modify list
            f5.close()

        #preprocessed file is at the tip so just read and give back
        else:
            #print(" tail is at the newest file")
            data = np.array(f5['processed'][:,:,:,:])
            f5.close()

            return data


    # actual file reading
    ret_arr = []
    for ff in files:
        arr = combine_tiles(ff, var, conf, isp=isp)
        
        if not(func == None):
            arr = func(arr)
        ret_arr.append( arr )
    data = np.array(ret_arr)


    # process quick read file

    #create it if it does not exist
    if not(os.path.isfile(fname_prepro)):
        #print(" creating preprocess file")
        f5 = h5py.File(fname_prepro,'w')
        dset = f5.create_dataset('processed', data=data)
        dset.attrs['current_head'] = last_file_num

        f5.close()

    #else update tail of processed file
    else:
        if newest_snapshot < last_file_num:
            #print("updating tail from {} to {}".format(newest_snapshot, last_file_num))
            f5 = h5py.File(fname_prepro,'r+')

            data_old = np.array(f5['processed'][:,:,:,:])
            #print(np.shape(data_old))
            #print(np.shape(data))

            data_new = np.append(data_old, data, axis=0)
            #print(np.shape(data_new))

            #update f5 file
            del f5['processed']
            dset = f5.create_dataset('processed', data=data_new)
            dset.attrs['current_head'] = last_file_num
            f5.close()

            return data_new


    return data



# combine tiles inside grid together into one array
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









