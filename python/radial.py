from distutils.fancy_getopt import fancy_getopt
import math
import numpy as np
import scipy as sp
import scipy.integrate
import scipy.sparse
import scipy.io
from scipy.sparse.linalg import (
    LinearOperator,
    gmres,
    lgmres,
    bicg,
    cgs,
    eigs,
    eigsh,
    spsolve,
    aslinearoperator,
    factorized,
)

from build.hamiltonian import PyHermiteBC


class RadialEquation:
    def __init__(self, grid_x, left_bc=1, right_bc=1):
        self.m1 = 1836.0
        self.m2 = 1.0
        self.format = "csc"

        self.grid_x = grid_x
        self.spl_x = PyHermiteBC(grid_x, left_bc, right_bc)

        self.ident = self.spl_x.getSPMatr()

        self.nx = self.spl_x.getSNMatr()
        self.b = self.ident

        self.nMatr = self.nx
        self.init_laplacian()
        self.init_potential()

    def init_laplacian(self):
        self.ddx = np.array(
            [
                [self.spl_x.d2BSplineBC(xi, i) for i in range(self.spl_x.getSpaceDim())]
                for j, xi in enumerate(self.spl_x.collocPoints)
            ]
        )

        self.ddx = sp.sparse.csc_matrix(self.ddx)
        m1, m2 = self.m1, self.m2
        mu_x = self.mu_x = 2.0 * m1 * m2 / (m1 + m2)
        print(f"reduced mass: {mu_x}")
        print("laplace is ready")

        self.laplace = 1.0 / mu_x * self.ddx

    def init_potential(self):
        self.lv = lambda x: 1.0 / x

        potential = np.array([self.lv(xi) for xi in self.spl_x.collocPoints])
        self.v_colloc = sp.sparse.spdiags(
            potential, 0, potential.size, potential.size, format=self.format
        )
        self.v = self.v_colloc @ self.ident
        print("potential is ready")

    def init_hamiltonian(self):
        s_inv = self.s_inv = sp.sparse.linalg.inv(self.ident)
        ddx_new = self.ddx_new = 1 / self.mu_x * self.ddx @ s_inv
        row = np.squeeze(ddx_new[0, :].toarray())
        col = np.squeeze(ddx_new[:, 0].toarray())
        # print(row.shape)
        # print(col.shape)

        diag = np.zeros(ddx_new.shape[0])
        diag[0] = diag[1] = 1
        for i in range(2, ddx_new.shape[0]):
            diag[i] = row[i] / col[i]
        p = np.diag(diag)
        p_inv = np.diag(1 / diag)
        self.p = p
        print(p)
        print(-p @ self.ddx_new)
        self.h = -p @ self.ddx_new - p @ self.v_colloc

    def init_default(self):
        self.h = -self.laplace - self.v

    def splinef(self, evec, x=1.0):
        values = np.zeros(self.spl_x.getSpaceDim())

        for j in range(self.spl_x.getSpaceDim()):
            if self.spl_x.fBSplineBC(x, j) == 0.0:
                continue
            values[j] = self.spl_x.fBSplineBC(x, j)

        return evec.dot(values)

    def eigs(self, **kwargs):
        evals, evecs = eigsh(self.h, M=self.p, **kwargs)

        self.eigenvalues, self.eigenvectors = evals, evecs

        return evals, evecs

    def exp_hamiltonian(self, dt=1e-2, **kwargs):
        inv_b = sp.sparse.linalg.inv(self.b)

        forward = self.b - 0.5 * dt * self.h
        backward = self.b + 0.5 * dt * self.h

        scipy.sparse.linalg.use_solver(useUmfpack=False)
        # mv = lambda x: backward @ x
        # op = LinearOperator(backward.shape, matvec=mv)

        solve = factorized(backward)
        mv = lambda x: solve(forward @ x)
        op = LinearOperator(backward.shape, matvec=mv)
        evals, evecs = eigs(op, **kwargs)
        evals = -np.log(evals) / dt
        for i in range(len(evals)):
            norm = sp.dot(evecs[:, i], self.nMatr @ evecs[:, i])
            evecs[:, i] = evecs[:, i] / math.sqrt(norm)
            phasenorm = self.splinef(evecs[:, i], 1.0)
            phasenorm = phasenorm / abs(phasenorm)
            evecs[:, i] = evecs[:, i] / phasenorm
        self.eigenvalues, self.eigenvectors = evals, evecs
        return evals, evecs
