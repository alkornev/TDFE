import build.hamiltonian
from build.hamiltonian import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import datetime
from matplotlib.animation import ArtistAnimation, FuncAnimation
from typing import Callable, Tuple
from tqdm import tqdm
import sys


BC = [-1, -1, -1, -1]
MASSES = [1836.0, 1836.0, 1.0]
REG = 0.25
LETTERS = [r"\alpha", r"\beta", r"\gamma"]
N_PLOTS = 2


def total_potential(x, y):
    potential = 1/np.sqrt(x*x+REG) - \
                1/np.sqrt((y-MASSES[0]*x/(MASSES[0]+MASSES[1]))**2+REG) - \
                1/np.sqrt((y+MASSES[1]*x/(MASSES[0]+MASSES[1]))**2+REG)
    return potential


def grid(x, R=2.0):
    if np.abs(x) < R:
        return R*np.sin(x/R*np.pi/2)
    else:
        return x


class AnimationWF:
    def __init__(self, box_shapes: Tuple[float], hamiltonian):
        aLeftEnd = box_shapes[0]
        aRightEnd = box_shapes[1]
        bLeftEnd = box_shapes[2]
        bRightEnd = box_shapes[3]

        self.xs = np.linspace(aLeftEnd, aRightEnd, 100)
        self.ys = np.linspace(bLeftEnd, bRightEnd, 100)
        self.X, self.Y = np.meshgrid(self.xs, self.ys)

        self.fig, self.ax = plt.subplots(1, 3, figsize=(15,7))
        self.fig.tight_layout(pad=3.0)
        self.text = self.fig.text(1.5, 0.04, 'step: 0', ha='center')

        self.cmap = 'RdBu_r'


        coefs = hamiltonian.get_eigenvectors()
        self.state = coefs[:, 0]
        self.h = hamiltonian

        self.dt = 0.001
        self.V = 0.075

    def update(self, i):
        self.state = self.h.evolutionStep(self.state, i, self.dt, self.V)            
        Z_real = np.array([[self.h.get_state(np.real(self.state), xi, yi) for xi in self.xs] for yi in self.ys])
        Z_imag = np.array([[self.h.get_state(np.imag(self.state), xi, yi) for xi in self.xs] for yi in self.ys])
        Z_ampl = np.array([[self.h.get_state(np.real(self.state), xi, yi)**2 + self.h.get_state(np.imag(self.state)**2, xi, yi) for xi in self.xs] for yi in self.ys])
        
        self.clear()
        self.text = self.fig.text(0.5, 0.04, f'step: {i}', ha='center')

        self.ax[0].contourf(self.X, self.Y, Z_real, vmin=-0.12, vmax=0.12, cmap=self.cmap)
        self.ax[1].contourf(self.X, self.Y, Z_imag, vmin=-0.12, vmax=0.12, cmap=self.cmap)
        self.ax[2].contourf(self.X, self.Y, Z_ampl, vmin=-0.12, vmax=0.12, cmap=self.cmap)

    def clear(self):
        for i in range(3):
            self.ax[i].clear()
        self.text.remove()
        self.ax[0].title.set_text('Real Part')
        self.ax[1].title.set_text('Imaginary Part')
        self.ax[2].title.set_text('Squared Amplitude')

    def animate(self, n_frames, interval):
        self.frames = n_frames
        #self.duration = self.dt*frames
        self.anim = FuncAnimation(self.fig, self.update, frames=n_frames, interval=interval, cache_frame_data=False, save_count=n_frames)
        self.save()

    def save(self):
        time = datetime.datetime.now().strftime("%m-%d-%Y-%H-%M-%S")
        self.anim.save(f'mp4s/WF-{time}-steps-{self.frames}.mp4', dpi=80)



def main():
    print("Module: ", dir(build.hamiltonian))
    set_num_threads(8)
    print("max threads: ", get_max_threads())
    aLeftEnd = -5
    aRightEnd = 5
    bLeftEnd = -5
    bRightEnd = 5
    box_shapes = (aLeftEnd, aRightEnd, bLeftEnd, bRightEnd)

    UnitAGrid = aRightEnd*np.linspace(-1, 1, 30)
    UnitBGrid = np.linspace(-1, 1, 30)

    aGrid = np.fromiter(map(lambda x: grid(x), UnitAGrid), dtype=np.float64)
    bGrid = bRightEnd*UnitBGrid


    begin_time = datetime.datetime.now()

    h3 = PyHamiltonian3D(aGrid, bGrid, BC, MASSES, REG, 0)
    print(h3.pMatr.shape)
    h3.get_spectrum()
    spectra = h3.get_eigenvalues()
    n = 0
    print(f"ground state: ", spectra[n])
    
    print("=="*40)

    anim = AnimationWF(box_shapes, h3)
    anim.animate(n_frames=2, interval=10)
        #cmap = 'RdBu_r
        #pair_indices = fr"{LETTERS[i%3]}"
        #triple_indices = fr"{LETTERS[(i+1)%3]}{LETTERS[(i+2)%3]}"

        #cb1 = ax[0, i].contourf(X, Y, Z_real, vmin=-0.12, vmax=0.12, cmap=cmap)
        #cb2 = ax[0, 1].contourf(X, Y, Z_imag, vmin=-0.12, vmax=0.12, cmap=cmap)

        #ax[0, i].set_xlabel(f'$r_{{{pair_indices}}}$', fontsize=20)
        #ax[0, i].set_ylabel(f'$r_{{{triple_indices}}}$', fontsize=20)
        #fig.colorbar(cb1, ax=ax[0, i])
        #fig.colorbar(cb2, ax=ax[0, 1])
        #ax[i].set_aspect('equal')
    print("Elapsed time:", datetime.datetime.now() - begin_time)
    #plt.show()
# def main():
#     print("Module: ", dir(build.hamiltonian))
#     set_num_threads(8)
#     print("max threads: ", get_max_threads())
#     aLeftEnd = -4
#     aRightEnd = 4
#     bLeftEnd = 10
#     bRightEnd = 10
#     box_shapes = (aLeftEnd, aRightEnd, bLeftEnd, bRightEnd)

#     UnitAGrid = aRightEnd*np.linspace(-1, 1, 20)
#     UnitBGrid = np.linspace(-1, 1, 10)**3

#     aGrid = np.fromiter(map(lambda x: grid(x), UnitAGrid), dtype=np.float64)
#     bGrid = bRightEnd*UnitBGrid


#     begin_time = datetime.datetime.now()

#     for i in range(1):
#         x_conv = []
#         y_conv = []
#         for size in np.linspace(15, 16, 1):
#             aRGrid = aRightEnd*np.linspace(-1, 1, int(size))
#             aGrid = np.fromiter(map(lambda x: grid(x), aRGrid), dtype=np.float64)
#             h3 = PyHamiltonian3D(aGrid, bGrid, BC, MASSES, REG, i)
#             h3.get_spectrum()
#             spectra = h3.get_eigenvalues()
#             coefs = h3.get_eigenvectors()
#             n = 0
#             print(f"{int(size)} ground state: ", spectra[n])
#             x_conv.append(1/int(size)**4)
#             y_conv.append(spectra[n])
#             dt = 0.001
#             coef = coefs[:, n]
            
#             print("=="*40)
#             #print(coef)

#         ax[1, i].scatter(x_conv, y_conv)
#         ax[1, i].set_xlim(0, 1e-4)

#         cmap = 'RdBu_r'
#         pair_indices = fr"{LETTERS[i%3]}"
#         triple_indices = fr"{LETTERS[(i+1)%3]}{LETTERS[(i+2)%3]}"

#         cb1 = ax[0, i].contourf(X, Y, Z_real, vmin=-0.12, vmax=0.12, cmap=cmap)
#         cb2 = ax[0, 1].contourf(X, Y, Z_imag, vmin=-0.12, vmax=0.12, cmap=cmap)

#         ax[0, i].set_xlabel(f'$r_{{{pair_indices}}}$', fontsize=20)
#         ax[0, i].set_ylabel(f'$r_{{{triple_indices}}}$', fontsize=20)
#         fig.colorbar(cb1, ax=ax[0, i])
#         fig.colorbar(cb2, ax=ax[0, 1])
#         #ax[i].set_aspect('equal')
#     print("Elapsed time:", datetime.datetime.now() - begin_time)
#     plt.show()

    

if __name__ == '__main__':
    main()
