import numpy as np

def sample_oscillating_langevin_antenna(size: int,
                                        characteristic_freq: float,
                                        decorrelation_rate: float,
                                        a0: complex = 1,
                                        dt: float = 1,
                                        gen = np.random.default_rng()):
    """
    Generates `size` amount of complex values
    according to oscillating Langevin antenna (TenBarge et al. 2014 eq. 10).

    Generated values can be used to drive antenna_mode time evolution,
    by giving them to lap_coeffs argument of antenna_mode.
    """
    if size == 0:
        return np.array([])

    a = np.zeros((size), dtype=complex)
    a[0] = a0

    w0, g0 = characteristic_freq, decorrelation_rate
    wa = w0 + 1j * g0

    u_re, u_im = gen.random(size - 1) - 1/2, gen.random(size - 1) - 1/2
    u = u_re + 1j * u_im

    sigma = np.abs(a0) * np.sqrt(12 * np.abs(g0) / dt)

    for n in range(1, size):
        F = sigma * u[n - 1] * dt
        a[n] = a[n - 1] * np.exp(-1j * wa * dt) + F

    return a


if __name__ == "__main__":
    def plot_oscillating_langevin_antenna_sample():
        import matplotlib.pyplot as plt

        def complex_scatter_plot(arr):
            re, im = np.real(arr), np.imag(arr)

            fix, ax = plt.subplots(1)
            ax.scatter(re, im, s=0.1, c=np.linspace(0, 1, len(re)))

            xlim = np.max(np.abs(re))
            ylim = np.max(np.abs(im))
            lim = np.max((xlim, ylim, 1))
            lims = (-lim, lim)

            ax.set(xlim=lims, ylim=lims, box_aspect=1)

            ax.axhline(y=0)
            ax.axvline(x=0)

            # Draw unit circle:
            theta = np.linspace(0, 2 * np.pi, 100)
            ax.plot(np.cos(theta), np.sin(theta))

            return fix, ax


        a = sample_oscillating_langevin_antenna(10000, 0.01, -0.0001)
        complex_scatter_plot(10 * a)

        plt.show()

    plot_oscillating_langevin_antenna_sample()
