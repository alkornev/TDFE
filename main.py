import numpy as np

from python.Hamiltonian1D3B import Hamiltonian1D3B
from python.plot import AnimationWF


gx = 10 * np.linspace(0, 1, 20) ** 2
gy = 50 * np.linspace(-1, 1, 20) ** 3


h = Hamiltonian1D3B(gx, gy)

evals, evecs = h.exp_hamiltonian()
h.evolutionStep(1, )

aLeftEnd = 0
aRightEnd = 10
bLeftEnd = -50
bRightEnd = 50
box_shapes = (aLeftEnd, aRightEnd, -25, 25)
plot = AnimationWF(box_shapes, h, 0.01, evecs[:, 0])
plot.plot(show=True)

