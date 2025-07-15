import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    for i in range(40):
        Ex = np.load(f"Ex{i}.npy")
        Ey = np.load(f"Ey{i}.npy")
        Ez = np.load(f"Ez{i}.npy")
        Bx = np.load(f"Bx{i}.npy")
        By = np.load(f"By{i}.npy")
        Bz = np.load(f"Bz{i}.npy")

        fig, (axE, axB) = plt.subplots(2, 1)

        axE.quiver(Ex[10:20, 10:20, 20], Ey[10:20, 10:20, 20], scale=8, label="E")
        axB.quiver(Bx[10:20, 10:20, 20], By[10:20, 10:20, 20], scale=8, label="B")

        axE.legend()
        axB.legend()

        plt.show()

