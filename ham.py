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


BC = [1, 1, 1, 1]
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
    def __init__(self, box_shapes: Tuple[float], hamiltonian, dt):
        aLeftEnd = box_shapes[0]
        aRightEnd = box_shapes[1]
        bLeftEnd = box_shapes[2]
        bRightEnd = box_shapes[3]

        self.xs = np.linspace(aLeftEnd, aRightEnd, 100)
        self.ys = np.linspace(bLeftEnd, bRightEnd, 100)
        self.X, self.Y = np.meshgrid(self.xs, self.ys)

        self.fig, self.ax = plt.subplots(figsize=(15,7))
        self.fig.tight_layout(pad=3.0)
        self.text = self.fig.text(1.5, 0.04, 'step: 0', ha='center')
        self.cb_ax = self.fig.add_axes([0.83, 0.1, 0.02, 0.8])
        self.cmap = 'RdBu_r'
        plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)


        coefs = hamiltonian.get_eigenvectors()
        self.state = coefs[:, 0]
        self.h = hamiltonian

        self.dt = dt
        self.V = 20.0

    def plot(self, show=False):        
        Z_ampl = np.array([[self.h.get_state(np.real(self.state), xi, yi)**2 + self.h.get_state(np.imag(self.state), xi, yi)**2 for xi in self.xs] for yi in self.ys])
        cf = self.ax.contourf(self.X, self.Y, Z_ampl, vmin=-0.01, vmax=1.5, cmap=self.cmap, levels=np.arange(-0.01, 1.52, 0.05))
        if show:
            plt.show()
        return cf

    def update(self, i):
        if i % 20 == 0:
            self.clear()
            self.text = self.fig.text(0.5, 0.04, f'step: {i}, time {self.dt*i:.2f}', ha='center')

            cf = self.plot()
            self.cbar = self.fig.colorbar(cf, cax=self.cb_ax)
        self.state = self.h.evolutionStep(self.state, i, self.dt, self.V)            

    def clear(self):
        for i in range(3):
            self.ax.clear()
        self.text.remove()
        self.ax.title.set_text('Squared Amplitude')

    def animate(self, n_frames, interval):
        self.frames = n_frames
        #self.duration = self.dt*frames
        self.anim = FuncAnimation(self.fig, self.update, frames=n_frames, interval=interval, cache_frame_data=False, save_count=100)
        self.save()

    def save(self):
        time = datetime.datetime.now().strftime("%m-%d-%Y-%H-%M-%S")
        self.anim.save(f'mp4s/WF-{time}-steps-{self.frames}.mp4', dpi=80)



def main():
    print("Module: ", dir(build.hamiltonian))
    set_num_threads(8)
    print("max threads: ", get_max_threads())
    aLeftEnd = 0
    aRightEnd = 10
    bLeftEnd = -45
    bRightEnd = 45
    box_shapes = (aLeftEnd, aRightEnd, bLeftEnd, bRightEnd)

    UnitAGrid = np.linspace(0, 1, 15)
    UnitBGrid = np.linspace(-1, 1, 45)

    aGrid = aRightEnd*UnitAGrid**3#aRightEnd*np.fromiter(map(lambda x: grid(x), UnitAGrid), dtype=np.float64)
    bGrid = bRightEnd*UnitBGrid


    begin_time = datetime.datetime.now()

    h3 = PyHamiltonian3D(aGrid, bGrid, BC, MASSES, REG, 0)
    dt = 0.01
    print(type(dt))
    h3.initLU(dt)
    print(h3.pMatr.shape)
    h3.get_spectrum(1, 50)
    spectra = h3.get_eigenvalues()
    n = 0
    print(f"ground state: ", spectra)
    
    print("=="*40)

    anim = AnimationWF(box_shapes, h3, dt)
    #anim.plot(show=True)
    anim.animate(n_frames=5000, interval=20)

    print("Elapsed time:", datetime.datetime.now() - begin_time)
    

if __name__ == '__main__':
    main()
