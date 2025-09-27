import numpy as np
from scipy.special import kn

def sample_juttner_synge_sobol_method(size: int,
                                      theta: float,
                                      gen = np.random.default_rng()):
    """
    Generates `size` amount of floats distributed
    according to Juttner-Synge (Zenitani 2015 eq. 2),
    with `theta = kT / (mc^2)` using sobol method.

    Random numbers are generated with the given generator.
    """

    eff = kn(2, 1 / theta) / (2 * theta**2)
    gen_coeff = 1.001 # Maybe could be adjusted for minor performance increase?

    candidates = []
    total_num_of_candidates = 0

    while total_num_of_candidates < size:
        N = int(gen_coeff * (size - total_num_of_candidates) / eff)
        x4 = gen.random(size=N)
        x5 = gen.random(size=N)
        x6 = gen.random(size=N)
        x7 = gen.random(size=N)

        u = -theta * np.log(x4*x5*x6)
        n = -theta * np.log(x4*x5*x6*x7)

        mask = (n**2 - u**2) > 1.0
        candidates.append(u[mask])
        total_num_of_candidates += len(candidates[-1])

    return np.concatenate(candidates)[:size]


def sample_maxwellian_box_muller_method(size: int,
                                        theta: float,
                                        gen = np.random.default_rng()):
    """
    Generates `size` amount of floats distributed
    according to Maxwell-Boltzmann distribution,
    with `theta = kT / m` using Box-Muller method (*).

    Random numbers are generated with the given generator.

    (*) Technically this is inverse transfrom sampling-method for Rayleigh distribution.
    """

    vth = np.sqrt(2 * theta)
    return vth * np.sqrt(-2.0 * np.log(gen.random(size=size)))


def sample_boosted_juttner_synge(size: int,
                                 theta: float,
                                 beta: None | float = None,
                                 Gamma: None | float = None,
                                 direction: str = '-x',
                                 gen = np.random.default_rng()):
    """
    Generates `size` amount of 3D spatial components of 4-velocities
    according to boosted Juttner-Synge distribution (Zenitani 2015),
    with `theta = kT / (mc^2)`.

    Boosted frame moves with 4-velocity (Gamma, -Gamma * beta).
    Either `Gamma` or `beta` can be specified but not both.

    `direction` defines direction of spatial part of the boost:
    '-x' => (-Gamma * beta, 0, 0)
     'x' => (Gamma * beta, 0, 0)
    '-y' => (0, -Gamma * beta, 0)
     'y' => (0, Gamma * beta, 0)
    '-z' => (0, 0, -Gamma * beta)
     'z' => (0, 0, Gamma * beta)

    Random numbers are generated with the given generator.

    Returns the components as a three 1D numpy arrays: (ux, uy, uz)
    """

    beta_given, Gamma_given = beta != None, Gamma != None

    if beta_given == Gamma_given:
        raise RuntimeError("Either beta or Gamma has to be given but not both.")

    if beta_given:
        Gamma = 1.0 / np.sqrt(1.0 - beta**2)
    else:
        beta = np.sqrt(1 - Gamma**(-2))

    if theta > 0.2:
        u = sample_juttner_synge_sobol_method(size, theta, gen)
    else:
        # Use Box-Muller for non-relativistic limit.
        u = sample_maxwellian_box_muller_method(size, theta, gen)

    X1 = gen.random(size=size)
    X2 = gen.random(size=size)
    ux = u * (2 * X1 - 1)
    uy = 2 * u * np.sqrt(X1 * (1 - X1)) * np.cos(2 * np.pi * X2)
    uz = 2 * u * np.sqrt(X1 * (1 - X1)) * np.sin(2 * np.pi * X2)

    vx = ux / np.sqrt(1.0 + u**2)
    flip_mask = (-beta * vx) > gen.random(size=size)
    ux[flip_mask] *= -1

    ux = Gamma * (ux + beta * np.sqrt(1 + u**2))

    match direction:
        case "-x":
            return ux, uy, uz
        case "x":
            return -ux, uy, uz
        case "-y":
            return uz, ux, uy
        case "y":
            return uz, -ux, uy
        case "-z":
            return uy, uz, ux
        case "z":
            return uy, uz, -ux
        case _:
            raise RuntimeError("Invalid direction given!")


if __name__ == "__main__":
    def plot_zenitani_2015_test_problem(res: int = 200,
                                        theta: float = 1,
                                        number_of_samples: int = int(1e6),
                                        Gammas = [1, 1.1, 10]):
        """
        Plots FIG. 3 from Zenitani 2015 using batch sampling
        and old sampling from pytools.
        """
        import time
        import matplotlib.pyplot as plt
        plt.rcParams['text.usetex'] = True
        from pytools import sample_boosted_maxwellian


        def expected_distribution(u, T: float, Gamma: float):
            """
            Eq. 33 of Zenitani 2015
            """
            a = Gamma * np.sqrt(1 + u**2) + T
            b = 2 * Gamma**2 * kn(2, 1 / T)

            beta = np.sqrt(1 - Gamma**(-2))
            c = Gamma * (np.sqrt(1 + u**2) - beta * u)
            return (a / b) * np.exp(-c / T) / Gamma


        u_arr = np.linspace(-20, 50, res)

        fig, (ax, ax_pytools) = plt.subplots(2, 1, layout="constrained")
        fig.suptitle(f"$\\theta = \\frac{{kT}}{{mc^2}} = {theta}, N = {number_of_samples}$")

        for G in Gammas:
            expected = expected_distribution(u_arr, theta, G)
            ax.plot(u_arr, expected, label=f"$\\Gamma={G}$")
            ax_pytools.plot(u_arr, expected, label=f"$\\Gamma={G}$")

            begin = time.time()
            ux, _, _, = sample_boosted_juttner_synge(number_of_samples, theta, Gamma=G)
            end = time.time()
            print(f"new: {end - begin:.3}s")

            begin = time.time()
            ux_pytools = [sample_boosted_maxwellian(theta, G, direction=1, dims=3)[0]
                          for _ in range(number_of_samples)]
            end = time.time()
            print(f"pytools: {end - begin:.3}s")

            ax.hist(ux, bins=100, density=True, alpha=0.2)
            ax_pytools.hist(ux_pytools, bins=100, density=True, alpha=0.2)

        ax.legend()
        ax.set(xlabel="$u'_x$",
               yscale="log",
               ylim=(1e-5, 1),
               xlim=(-20, 50),
               title="new")

        ax_pytools.legend()
        ax_pytools.set(xlabel="$u'_x$",
                       yscale="log",
                       ylim=(1e-5, 1),
                       xlim=(-20, 50),
                       title="pytools")
        plt.show()

    plot_zenitani_2015_test_problem()
