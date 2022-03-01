import datetime
import sys
from typing import Callable, Tuple

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.sparse.linalg
from matplotlib import cm
from matplotlib.animation import ArtistAnimation, FuncAnimation
from tqdm import tqdm


class AnimationWF:
    def __init__(self, box_shapes: Tuple[float], hamiltonian, dt, coefs):
        aLeftEnd = box_shapes[0]
        aRightEnd = box_shapes[1]
        bLeftEnd = box_shapes[2]
        bRightEnd = box_shapes[3]

        self.xs = np.linspace(aLeftEnd, aRightEnd, 100)
        self.ys = np.linspace(bLeftEnd, bRightEnd, 100)
        self.X, self.Y = np.meshgrid(self.xs, self.ys)

        self.fig, self.ax = plt.subplots(figsize=(15, 7))
        self.fig.tight_layout(pad=3.0)
        self.text = self.fig.text(1.5, 0.04, "step: 0", ha="center")
        self.cb_ax = self.fig.add_axes([0.83, 0.1, 0.02, 0.8])
        self.cmap = "RdBu_r"
        plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)

        self.state = coefs
        self.h = hamiltonian

        self.dt = dt

    def plot(self, show=False):
        Z_ampl = np.array(
            [
                [
                    self.h.splinef(np.real(self.state), xi, yi) ** 2
                    + self.h.splinef(np.imag(self.state), xi, yi) ** 2
                    for xi in self.xs
                ]
                for yi in self.ys
            ]
        )
        cf = self.ax.contourf(
            self.X,
            self.Y,
            Z_ampl,
            vmin=-0.01,
            vmax=1.3,
            cmap=self.cmap,
            levels=np.arange(-0.01, 1.3, 0.1),
            extend="max"
        )
        if show:
            plt.show()
        return cf
    
    def update(self, i):
        if i % 5 == 0:
            self.clear()
            self.text = self.fig.text(
                0.5, 0.04, f"step: {i}, time {self.dt*i:.2f} fs FWHM", ha="center"
            )

            cf = self.plot()
            self.cbar = self.fig.colorbar(cf, cax=self.cb_ax)
        self.state = self.h.evolutionStep(self.state, i)

    def clear(self):
        for i in range(3):
            self.ax.clear()
        self.text.remove()
        self.ax.title.set_text("Squared Amplitude")

    def animate(self, n_frames, interval):
        self.frames = n_frames
        # self.duration = self.dt*frames
        self.anim = FuncAnimation(
            self.fig,
            self.update,
            frames=n_frames,
            interval=interval,
            cache_frame_data=False,
            save_count=100,
        )
        self.save()

    def save(self):
        time = datetime.datetime.now().strftime("%m-%d-%Y-%H-%M-%S")
        self.anim.save(f"mp4s/WF-{time}-steps-{self.frames}.mp4", dpi=80)
