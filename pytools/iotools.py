import numpy as np


# read simulation output file and reshape to python format
def read_h5_array(f5, var_name, stride=1):

    try:
        nx = f5['Nx'][()]
        ny = f5['Ny'][()]
        nz = f5['Nz'][()]
    except:
        nx,ny,nz = np.shape(f5[var_name])
        print("read_h5_array: fallback to brute-force reading...")
        val = f5[var_name][()]
        return val

    # column-ordered data with image convention (x horizontal, y vertical, ..)
    # i.e., so-called fortran ordering
    if stride == 1:
        val = f5[var_name][()]
    else:
        val = f5[var_name][::stride**3]
        nx = int(nx/stride)
        ny = int(ny/stride)
        nz = int(nz/stride)

    #print("reshaping 1D array of {} into multiD with {} {} {}".format(len(val), nx,ny,nz))

    # reshape to python format; from Fortran image to C matrix
    val = np.reshape(val, (nz, ny, nx))
    val = val.ravel(order='F').reshape((nx,ny,nz))

    return val


