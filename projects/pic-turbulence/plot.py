import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


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

def read_full_box(path_to_h5: str, var: str):
    f5 = h5py.File(path_to_h5, "r")
    return read_h5_array(f5, var)


def read_je(path_to_h5: str):
    jx = read_full_box(path_to_h5, "jx")
    jy = read_full_box(path_to_h5, "jy")
    jz = read_full_box(path_to_h5, "jz")
    ex = read_full_box(path_to_h5, "ex")
    ey = read_full_box(path_to_h5, "ey")
    ez = read_full_box(path_to_h5, "ez")

    return jx * ex + jy * ey + jz * ez


def plot(data, fig, ax_xy, ax_yz, ax_zx, vmin=None, vmax=None, cblabel=""):

    x_plane, y_plane, z_plane = 160, 20, 3 
    xy_data = data[:, :, z_plane]
    yz_data = data[x_plane, :, :]
    zx_data = data[:, y_plane, :]
    vmin = vmin if vmin else min(np.min(xy_data), np.min(yz_data), np.min(zx_data))
    vmax = vmax if vmax else max(np.max(xy_data), np.max(yz_data), np.max(zx_data))
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    ax_xy.imshow(xy_data, norm=norm)
    ax_xy.set_title(f"z={z_plane}")

    ax_yz.imshow(yz_data, norm=norm)
    ax_yz.set_title(f"x={x_plane}")

    img = ax_zx.imshow(zx_data, norm=norm)
    ax_zx.set_title(f"y={y_plane}")

    fig.colorbar(img,
                 ax=[ax_xy, ax_yz, ax_zx],
                 orientation="vertical",
                 label=cblabel)

if __name__ == "__main__":

    var = sys.argv[1]

    rows = len(sys.argv[2:])
    fig, ax_all = plt.subplots(rows, 3,
                               layout="compressed",
                               figsize=(2 * 3, 2 * rows))

    if rows == 1:
        ax_all = [ax_all]

    for i, file in enumerate(sys.argv[2:]):
        ax = ax_all[i]

        match var:
            case "je":
                plot(read_je(file), fig, *ax, cblabel=file)
            case "bx":
                plot(read_full_box(file, "bx"), fig, *ax, cblabel=file)
            case "by":
                plot(read_full_box(file, "by"), fig, *ax, cblabel=file)
            case "bz":
                plot(read_full_box(file, "bz"), fig, *ax, cblabel=file)
            case "ex":
                plot(read_full_box(file, "ex"), fig, *ax, cblabel=file)
            case "ey":
                plot(read_full_box(file, "ey"), fig, *ax, cblabel=file)
            case "ez":
                plot(read_full_box(file, "ez"), fig, *ax, cblabel=file)
            case "jx":
                plot(read_full_box(file, "jx"), fig, *ax, cblabel=file)
            case "jy":
                plot(read_full_box(file, "jy"), fig, *ax, cblabel=file)
            case "jz":
                plot(read_full_box(file, "jz"), fig, *ax, cblabel=file)

    fig.suptitle(var)
    plt.show()
