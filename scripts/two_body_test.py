import build.hamiltonian
from build.hamiltonian import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import datetime


LEFT_END = 0
RIGHT_END = 40
TWO_BODY_LEFT_BC = 0
TWO_BODY_RIGHT_BC = -1
TWO_BODY_MASSES = [1836.0, 1.0]
REG = 0.25

TO_PLOT_BASIS = True
TO_PLOT_BOUND_STATE = True
TO_PLOT_CONV = True


def main():
    print(dir(build.hamiltonian))
    print("max threads: ", get_max_threads())
    set_num_threads(8)

    x_conv = []
    y_conv = []

    fig, ax = plt.subplots(1, 3, figsize=(15, 7))
    fig.tight_layout(pad=3.0)

    print("Calculating the two-body ground state: ")
    for n in [60, 120, 180, 200, 300, 400]:
        h2_grid = RIGHT_END * np.linspace(-1, 1, n)
        splines = PyHermiteBC(h2_grid, TWO_BODY_LEFT_BC, TWO_BODY_RIGHT_BC)

        h2 = PyHamiltonian2D(
            h2_grid, TWO_BODY_LEFT_BC, TWO_BODY_RIGHT_BC, TWO_BODY_MASSES, REG
        )
        h2.get_spectrum()
        spectra = h2.get_eigenvalues()
        coefs = h2.get_eigenvectors()
        print(spectra[:2])

        x_conv.append(1 / n**4)

        y_conv.append(spectra[0])

    if TO_PLOT_BASIS:
        grd = (RIGHT_END + LEFT_END) * np.linspace(-1, 1, 1000) - LEFT_END
        for i in range(splines.getSpaceDim()):
            x = [splines.fBSplineBC(t, i) for t in grd]
            ax[0].plot(grd, x)

    if TO_PLOT_BOUND_STATE:
        grd = (RIGHT_END + LEFT_END) * np.linspace(-1, 1, 1000) - LEFT_END
        y = [h2.get_state(coefs[:, 0], xi) for xi in grd]
        ax[1].plot(grd, y)

    if TO_PLOT_CONV:
        ax[2].scatter(x_conv, y_conv)
        ax[2].set_xlim(-2e-9, 2e-7)
        # ax[2].set_yscale('symlog')
        # ax[2].set_xscale('log')
    plt.show()


if __name__ == "__main__":
    main()
