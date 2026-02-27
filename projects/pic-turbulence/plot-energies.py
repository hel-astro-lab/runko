import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

if __name__ == "__main__":

    outdir = sys.argv[1]
    kinetic_energy_path = f"{outdir}/average_kinetic_energy.txt"
    B_energy_path = f"{outdir}/average_B_energy_density.txt"
    E_energy_path = f"{outdir}/average_E_energy_density.txt"

    fig, (axP, axB, axE) = plt.subplots(3, 1)

    pdata = np.loadtxt(kinetic_energy_path)
    Bdata = np.loadtxt(B_energy_path)
    Edata = np.loadtxt(E_energy_path)

    axP.plot(pdata[:, 0], pdata[:, 1], label="p0")
    axP.plot(pdata[:, 0], pdata[:, 2], label="p1")
    axB.plot(Bdata[:, 0], Bdata[:, 1])
    axE.plot(Edata[:, 0], Edata[:, 1])

    axP.legend()

    plt.show()
