import build.hamiltonian
from build.hamiltonian import *
import numpy as np
import scipy as sp
import scipy.sparse.linalg
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
    def __init__(self, box_shapes: Tuple[float], hamiltonian, dt, coefs):
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

        self.state = coefs
        self.h = hamiltonian

        self.dt = dt

    def plot(self, show=False):        
        Z_ampl = np.array([[self.h.get_state(np.real(self.state), xi, yi)**2 + self.h.get_state(np.imag(self.state), xi, yi)**2 for xi in self.xs] for yi in self.ys])
        cf = self.ax.contourf(self.X, self.Y, Z_ampl, vmin=-0.01, vmax=1., cmap=self.cmap, levels=np.arange(-0.01, 1.3, 0.05))
        if show:
            plt.show()
        return cf

    def update(self, i):
        if i % 50 == 0:
            self.clear()
            self.text = self.fig.text(0.5, 0.04, f'step: {i}, time {self.dt*i:.2f}', ha='center')

            cf = self.plot()
            self.cbar = self.fig.colorbar(cf, cax=self.cb_ax)
        self.state = self.h.evolutionStep(self.state, i)            

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
    bLeftEnd = -50
    bRightEnd = 50
    box_shapes = (aLeftEnd, aRightEnd, -50, 50)

    UnitAGrid = np.linspace(0, 1, 120)
    UnitBGrid = np.linspace(-1, 1, 120)

    
    aGrid = aRightEnd*UnitAGrid#**2
    #aRightEnd*np.fromiter(mapambda x: grid(x), UnitAGrid), dtype=np.float64)
    bGrid = bRightEnd*UnitBGrid

    begin_time = datetime.datetime.now()

    h3 = PyHamiltonian3D(aGrid, bGrid, BC, MASSES, REG, 0)

    #void initImpulse(double init_time, double init_phase, double intensity, double freq, double duration, double step);
    t0 = -450.0
    phase = 0.0
    I0 = 0.05
    w0 = 0.058
    tau = 248.0
    dt = 0.5

    h3.init_impulse(t0, phase, I0, w0, tau, dt)
    h3.init_absorption(200)
    h3.init_LU()

    print(h3.nMatr.shape)
    print(type(h3.h))
    #vals, vecs = sp.sparse.linalg.eigs(h3.pMatr, 10, h3 which="SI", tol=1E-12)
    #vals, vecs = sp.sparse.linalg.eig
    h3.get_spectrum(3, 100)
    spectra = h3.get_eigenvalues()
    vecs = h3.get_eigenvectors()
    print("spectra: ", spectra)
    #n = 0
    print("=="*40)
    # coefs = vecs[:, 1]
    # norm = np.sqrt(np.dot(h3.nMatr.dot(coefs), coefs))

    # coefs = coefs/norm
    # phase = h3.get_state(coefs, 1.0, 1.0)
    # coefs = coefs * phase/np.abs(phase)

    anim = AnimationWF(box_shapes, h3, dt, vecs[:, 0])
    #anim.plot(show=True)
    anim.animate(n_frames=0, interval=10)

    print("Elapsed time:", datetime.datetime.now() - begin_time)
    

if __name__ == '__main__':
    main()
