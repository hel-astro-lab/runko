import sys
import h5py
import pytools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


def read_full_box(path_to_h5: str, var: str):
    f5 = h5py.File(path_to_h5, "r")
    return pytools.read_h5_array(f5, var)


def read_je(path_to_h5: str):
    jx = read_full_box(path_to_h5, "jx")
    jy = read_full_box(path_to_h5, "jy")
    jz = read_full_box(path_to_h5, "jz")
    ex = read_full_box(path_to_h5, "ex")
    ey = read_full_box(path_to_h5, "ey")
    ez = read_full_box(path_to_h5, "ez")

    return jx * ex + jy * ey + jz * ez


def plot(data, path_to_h5, vmin=None, vmax=None):

    fig, (ax_xy, ax_yz, ax_zx) = plt.subplots(1, 3)

    x_plane, y_plane, z_plane = 30, 50, 60
    xy_data = data[:, :, z_plane]
    yz_data = data[x_plane, :, :]
    zx_data = data[:, y_plane, :]
    vmin = vmin if vmin else min(np.min(xy_data), np.min(yz_data), np.min(zx_data))
    vmax = vmax if vmax else max(np.max(xy_data), np.max(yz_data), np.max(zx_data))
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    img = ax_xy.imshow(xy_data, norm=norm)
    ax_xy.set_title(f"z={z_plane}")

    ax_yz.imshow(yz_data, norm=norm)
    ax_yz.set_title(f"x={x_plane}")

    ax_zx.imshow(zx_data, norm=norm)
    ax_zx.set_title(f"y={y_plane}")

    fig.colorbar(img, ax=[ax_xy, ax_yz, ax_zx], orientation="horizontal")
    fig.suptitle(f"{path_to_h5}")

    plt.show()

if __name__ == "__main__":

    var = sys.argv[1]

    for file in sys.argv[2:]:
        match var:
            case "je":
                plot(read_je(file), f"je: file", vmin=-0.5, vmax=0.5)
            case "bx":
                plot(read_full_box(file, "bx"), f"bx: {file}")
            case "by":
                plot(read_full_box(file, "by"), f"by: {file}")
            case "bz":
                plot(read_full_box(file, "bz"), f"bz: {file}")
            case "ex":
                plot(read_full_box(file, "ex"), f"ex: {file}")
            case "ey":
                plot(read_full_box(file, "ey"), f"ey: {file}")
            case "ez":
                plot(read_full_box(file, "ez"), f"ez: {file}")
