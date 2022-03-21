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
    spsolve,
    aslinearoperator,
    factorized,
)

from build.hamiltonian import PyHermiteBC


class Hamiltonian1D3B:
    def __init__(self, grid_x, grid_y, left_bc=1, right_bc=1):
        self.m1 = 1836.0
        self.m2 = 1836.0
        self.m3 = 1.0
        self.format = 'csc'

        self.grid_x = grid_x
        self.spl_x = PyHermiteBC(grid_x, left_bc, right_bc)
        self.spl_y = PyHermiteBC(grid_y, left_bc, right_bc)

        self.ident_x = self.spl_x.getSPMatr()
        self.ident_y = self.spl_y.getSPMatr()

        self.nx = self.spl_x.getSNMatr()
        self.ny = self.spl_y.getSNMatr()

        self.my_sx = np.array(
            [
                [self.spl_x.fBSplineBC(xi, i) for i in range(self.spl_x.getSpaceDim())]
                for j, xi in enumerate(self.spl_x.collocPoints)
            ]
        )
        self.my_sx = sp.sparse.csc_matrix(self.my_sx)
        self.Rx = 1.0

        self.nMatr = sp.sparse.kron(self.nx, self.ny, format=self.format)
        self.ident = sp.sparse.kron(self.ident_x, self.ident_y, format=self.format)
        self.init_laplacian()
        self.init_potential()
        self.init_hamiltonian()
        self.init_absorption()

    def init_laplacian(self):
        self.dx = np.array(
            [
                [self.spl_x.d2BSplineBC(xi, i) for i in range(self.spl_x.getSpaceDim())]
                for j, xi in enumerate(self.spl_x.collocPoints)
            ]
        )
        self.dy = np.array(
            [
                [self.spl_y.d2BSplineBC(yi, i) for i in range(self.spl_y.getSpaceDim())]
                for j, yi in enumerate(self.spl_y.collocPoints)
            ]
        )
        self.dx = sp.sparse.csc_matrix(self.dx)
        self.dy = sp.sparse.csc_matrix(self.dy)
        m1, m2, m3 = self.m1, self.m2, self.m3
        mu_x = self.mu_x = 2.0 * m1 * m2 / (m1 + m2)
        mu_y = self.mu_y = 2.0 * ((m1 + m2) * m3) / (m1 + m2 + m3)

        self.laplace_x = -1.0 / mu_x * sp.sparse.kron(self.dx, self.ident_y, format=self.format)
        self.laplace_y = -1.0 / mu_y * sp.sparse.kron(self.ident_x, self.dy, format=self.format)
        self.laplace = self.laplace_x + self.laplace_y

        inertio = np.array(
            [
                xi * xi
                for xi in self.spl_x.collocPoints
                for yi in self.spl_y.collocPoints
            ]
        )
        self.inertio = sp.sparse.spdiags(inertio, 0, inertio.size, inertio.size, format=self.format)

    def init_potential(self, reg_x=0.03, reg_y=0.25):
        self.reg_x = reg_x
        self.reg_y = reg_y
        m1, m2, m3 = self.m1, self.m2, self.m3
        self.lv = lambda x, y: 1.0 / math.sqrt(x * x + self.reg_x) - (
            1.0
            / math.sqrt(
                (y - x / (1.0 + m2 / m1)) * (y - x / (1.0 + m2 / m1)) + self.reg_y
            )
            + 1.0
            / math.sqrt(
                (y + x / (1.0 + m2 / m1)) * (y + x / (1.0 + m2 / m1)) + self.reg_y
            )
        )

        potential = np.array(
            [
                self.lv(xi, yi)
                for xi in self.spl_x.collocPoints
                for yi in self.spl_y.collocPoints
            ]
        )
        self.v = sp.sparse.spdiags(potential, 0, potential.size, potential.size, format=self.format)
        self.v = 0.5 * self.v @ self.ident

    def init_hamiltonian(self):
        self.h = self.laplace + self.v

    def gauss(self, t):
        vt = lambda x: self.I0*math.exp(-math.pow(t/self.tau, 2)) * math.cos(t * self.w0 + self.phase)
        diag_vt = np.array(
            [
                -vt(t)*yi
                for _ in self.spl_x.collocPoints
                for yi in self.spl_y.collocPoints
            ]
        )
        diag_vt = sp.sparse.spdiags(diag_vt, 0, diag_vt.size, diag_vt.size, format=self.format)

        return diag_vt @ self.ident


    def splinef(self, evec, x=1.0, y=0.0):
        values = np.zeros(self.spl_x.getSpaceDim()*self.spl_y.getSpaceDim())
        for i in range(self.spl_x.getSpaceDim()):
            if self.spl_x.fBSplineBC(x, i) == 0.0:
                continue
            for j in range(self.spl_y.getSpaceDim()):
                if self.spl_y.fBSplineBC(y, j) == 0.0:
                    continue
                values[i*self.spl_y.getSpaceDim() + j] = self.spl_x.fBSplineBC(x, i) * self.spl_y.fBSplineBC(y, j)

        return evec.dot(values)

    def init_absorption(self, y0=48, x0=14, kx=1.0, ky=1.0):
        absorb_y= lambda y: 0.0 if abs(y) < y0 else -1.0j * (ky * math.pow(y - y0, 2))
        absorb_x= lambda x: 0.0 if x < x0 else -1.0j * (kx *  math.pow(x - x0, 2))
        absorption = np.array(
            [
                absorb_x(xi) + absorb_y(yi)
                for xi in self.spl_x.collocPoints
                for yi in self.spl_y.collocPoints
            ]
        )
        self.absorb = sp.sparse.spdiags(absorption, 0, absorption.size, absorption.size, format=self.format)
        self.absorb = self.absorb @ self.ident

    def scale(self, t, t0, v=0.01, n_photon=10, **kwargs):
        v = np.sqrt(2 * n_photon * self.w0 / (self.mu_x / 2)) / np.max(self.grid_x)
        if t >= t0:
            Rx2 = 1.0 + math.pow(v * (t - t0), 4)
            Rx = math.sqrt(math.sqrt(Rx2))
            self.Rx = Rx
            d2Rx = 3 * math.pow(v, 4) * math.pow(t - t0, 2) / math.pow(Rx, 7)

            self.laplace = 1.0 / (Rx * Rx) * self.laplace_x
            self.laplace += self.laplace_y
            self.laplace += 0.25 * self.mu_x * Rx * d2Rx * self.inertio @ self.ident

            potential = np.array(
                [
                    self.lv(Rx * xi, yi)
                    for xi in self.spl_x.collocPoints
                    for yi in self.spl_y.collocPoints
                ]
            )

            self.v = sp.sparse.spdiags(potential, 0, potential.size, potential.size, format=self.format)
            self.v = 0.5 * self.v @ self.ident
            self.h = self.laplace + self.v

    def exp_hamiltonian(self, dt=1e-2, **kwargs):
        forward = self.ident - 0.5 * dt * self.h
        backward = self.ident + 0.5 * dt * self.h

        scipy.sparse.linalg.use_solver(useUmfpack=False)
        # mv = lambda x: backward @ x
        # op = LinearOperator(backward.shape, matvec=mv)
        print("init")
        
        solve = factorized(backward)
        print('lu')
        mv = lambda x: solve(forward @ x)
        op = LinearOperator(backward.shape, matvec=mv)
        print(type(backward))
        evals, evecs = eigs(op, **kwargs)
        evals = -np.log(evals) / dt
        print('log')
        for i in range(len(evals)):
            print(i)
            norm = sp.dot(evecs[:, i], self.nMatr @ evecs[:, i])
            evecs[:, i] = evecs[:, i] / math.sqrt(norm)
            phasenorm = self.splinef(evecs[:, i], 1.0, 0.0)
            phasenorm = phasenorm / abs(phasenorm)
            evecs[:, i] = evecs[:, i] / phasenorm
        self.eigenvalues, self.eigenvectors = evals, evecs
        return evals, evecs

    def init_gauss(self, I0=0.05, duration=10, w0=0.058, phase=0.0):
        self.I0 = I0/3.5095 # convert intensity to I0 * 10^-16
        self.tau = duration/(2.41889e-2*1.665109222) # convert femtoseconds to atomic units FWHM
        self.phase = phase
        self.w0 = w0

    def time_to_au(self, t):
        return t / (2.41889e-2*1.665109222)

    def evolution_step(self, evec, init_time = 0, dt=1e-3, iters=100, debug=False, **kwargs):
        scipy.sparse.linalg.use_solver(useUmfpack=False)

        cache = []
        for i in range(iters):
            if debug:
                print(f"iter: {i}, time {i*0.5}, {np.linalg.norm(evec.conjugate().dot(self.nMatr @ evec))}")
            if 'cache_step' in kwargs and i % kwargs['cache_step'] == 0:
                cache.append(evec)
            t = i*dt + 0.5*dt + init_time
            if t >= kwargs['t0']:
                self.scale(t, **kwargs)


            imp = lambda x: 1.0 if abs(t) < self.tau / 2 else 0.0
            p_forward = self.ident - 0.5j * dt * (self.gauss(t) * imp(t) + self.h + self.absorb)
            p_backward = self.ident + 0.5j * dt * (self.gauss(t) * imp(t) + self.h + self.absorb)
            p_factorized = factorized(p_backward)
            p_mv = lambda x: p_factorized(p_forward @ x)
            evec = p_mv(evec)

            

        return evec, cache

    def old_evol(self, evec, dt=1e-3, t0_scale=0.0, iters=100, debug=False):
        cache = []
        t = 0
        for i in range(iters):
            if t >= t0_scale:
                self.scale(t + dt * 0.5, t0_scale, 0.0025)
            forward = self.ident - 0.5j * dt * (self.h + self.absorb)
            backward = self.ident + 0.5j * dt * (self.h + self.absorb)

            mv = lambda x: spsolve(backward, forward @ x)
            op = LinearOperator(backward.shape, matvec=mv)

            if debug:
                print(np.linalg.norm(evec.conjugate().dot(self.nMatr @ evec)))
            evec = op(evec)

            cache.append(evec.copy())

        return evec, cache
