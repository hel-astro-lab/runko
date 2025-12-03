import matplotlib.pyplot as plt
import numpy as np
import itertools

ax = plt.figure().add_subplot(projection='3d')

Nx, Ny, Nz = 4, 5, 6
h = 2

def offset(i, j, k):
    whole_i_slices = min(max(0, i - h), Nx - 2 * h)
    skips = whole_i_slices * (Ny - 2 * h) * (Nz - 2 * h)

    if i in range(h, Nx - h):
        whole_j_slices = min(max(0, j - h), Ny - 2 * h)
        skips += whole_j_slices * (Nz - 2 * h)

        if k >= Nz - h and j in range(h, Ny - h):
            skips += Nz - 2 * h

    return i * (Ny * Nz) + j * Nz + k - skips

offsets = []

for i, j, k in itertools.product(range(Nx), range(Ny), range(Nz)):
    A_i = not (i < h)
    B_i = not (i >= Nx - h)

    A_j = not (j < h)
    B_j = not (j >= Ny - h)

    A_k = not (k < h)
    B_k = not (k >= Nz - h)

    in_middle =  A_i and B_i and A_j and B_j and A_k and B_k

    if in_middle:
        continue
    o = offset(i, j, k)
    offsets.append(o)
    ax.text(i, j, k, f"{o}")

offsets.sort()

for i in range(len(offsets) - 1):
    if offsets[i + 1] - offsets[i] != 1:
        print(i, i + 1, offsets[i], offsets[i + 1])


ax.set(xlim=(-1, Nx),
       ylim=(-1, Ny),
       zlim=(-1, Nz))

plt.show()
