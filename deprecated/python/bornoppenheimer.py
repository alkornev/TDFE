import datetime
import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

import build.hamiltonian
from build.hamiltonian import *

LEFT_END = 0
RIGHT_END = 60
TWO_BODY_LEFT_BC = 1
TWO_BODY_RIGHT_BC = 1
BODY_MASSES = [1836.0, 1836.0, 1.0]
REG = 0.25


def main():
    print(dir(build.hamiltonian))
    print("max threads: ", get_max_threads())
    set_num_threads(8)

    print("Calculating Born-Oppenheimer electronic state: ")
    Rs = np.arange(0.00, 10.0, 0.25)
    U = np.zeros((3, len(Rs)))
    for k, R in enumerate(Rs):
        for n in [200]:
            h2_grid = RIGHT_END * np.linspace(-1, 1, n)

            h2 = PyBornOppenHeimer2D(
                h2_grid, TWO_BODY_LEFT_BC, TWO_BODY_RIGHT_BC, BODY_MASSES, REG, R
            )
            h2.get_spectrum()

            spectra = h2.get_eigenvalues()
            coefs = h2.get_eigenvectors()
            print(
                f"R: {R:.2f},    grid size: {len(h2_grid)},    U: {spectra[:3] + 0.5/math.sqrt(R*R + REG)}"
            )
            U[0, k] = spectra[0] + 0.5 / math.sqrt(R * R + REG)
            U[1, k] = spectra[1] + 0.5 / math.sqrt(R * R + REG)
            U[2, k] = spectra[2] + 0.5 / math.sqrt(R * R + REG)
            grd = RIGHT_END * np.linspace(-1, 1, 500)
            y = [h2.get_state(coefs[:, 0], xi) for xi in grd]
            # plt.plot(grd, y)
            # plt.show()
    plt.plot(Rs, U[0, :], Rs, U[1, :], Rs, U[2, :])
    plt.show()


if __name__ == "__main__":
    main()
