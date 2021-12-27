import build.hamiltonian
from build.hamiltonian import *
import numpy as np


np.set_printoptions(suppress=True, precision=3)

grid = np.sin(np.linspace(0, 1, 4)*np.pi/2)
h = PyHamiltonian2D(grid, 1, 0, [1,1], 1)


t = np.array([[],[],[],[],[]])
pInv = np.linalg.inv(h.pMatr)
d2p = h.d2Matr @ pInv
p2p = h.potential @ pInv
print(h.pMatr)
print(d2p)
print(np.sum(d2p, axis=0))
print(p2p)