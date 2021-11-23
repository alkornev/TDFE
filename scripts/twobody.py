import build.hamiltonian
from build.hamiltonian import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import datetime
from mpl_toolkits.mplot3d import Axes3D


BC = [1, 1, 1, 1]
MASSES = [1836.0, 1836.0, 1.0]
REG = 1.0
LETTERS = [r"\alpha", r"\beta", r"\gamma"]
N_PLOTS = 3

def total_potential(x, y):
    potential = 1/np.sqrt(x*x+REG) - \
                1/np.sqrt((y-MASSES[0]*x/(MASSES[0]+MASSES[1]))**2+REG) - \
                1/np.sqrt((y+MASSES[1]*x/(MASSES[0]+MASSES[1]))**2+REG)
    return potential

def main():
    aLeftEnd = -5
    aRightEnd = 5

    bLeftEnd = -5
    bRightEnd = 5

    # draw potential surface

    aRGrid = np.linspace(-1, 1, 5)
    aGrid = aRGrid*aRightEnd
    bRGrid = np.linspace(-1, 1, 5)
    bGrid = bRGrid*bRightEnd

    xs = np.linspace(aLeftEnd, aRightEnd, 200)
    ys = np.linspace(bLeftEnd, bRightEnd, 200)
    X, Y = np.meshgrid(xs, ys)


    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    fig.tight_layout(pad=3.0)


    begin_time = datetime.datetime.now()
    h3 = PyHamiltonian3D(aGrid, bGrid, BC, MASSES, REG, 0)
    length = h3.pMatr.shape[0]
    s = 25
    for i in range(s, s+1, 2):
        coef = np.array([1 if x == i else 0 for x in range(length)])
        Z = np.array([[h3.get_state(coef, xi, yi) for xi in xs] for yi in ys])

        cmap = 'RdBu_r'
        pair_indices = fr"{LETTERS[i%3]}"
        triple_indices = fr"{LETTERS[(i+1)%3]}{LETTERS[(i+2)%3]}"
        cb = ax.plot_surface(X, Y, Z, cmap=cmap)
    ax.set_xlabel(f'$r_{{{pair_indices}}}$', fontsize=20)
    ax.set_ylabel(f'$r_{{{triple_indices}}}$', fontsize=20)
    plt.show()

    

if __name__ == '__main__':
    main()
