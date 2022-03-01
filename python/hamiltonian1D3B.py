from distutils.fancy_getopt import fancy_getopt
import math
import numpy as np
from pandas import factorize
import scipy as sp
import scipy.integrate
import scipy.sparse
import scipy.io
from scipy.sparse.linalg import LinearOperator, gmres, lgmres, bicg, cgs, eigs, spsolve, aslinearoperator, factorized

from build.hamiltonian import PyHermiteBC


class Hamiltonian1D3B:
    def __init__(self, grid_x, grid_y, left_bc=1, right_bc=1):
        self.m1 = 1836.0
        self.m2 = 1836.0
        self.m3 = 1.0

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

        self.nMatr = sp.sparse.kron(self.nx, self.ny)
        self.eye = sp.sparse.kron(self.ident_x, self.ident_y)
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
        
        self.laplace_x = -1.0 / mu_x * sp.sparse.kron(self.dx, self.ident_y)
        self.laplace_y = -1.0 / mu_y * sp.sparse.kron(self.ident_x, self.dy)
        self.laplace = self.laplace_x + self.laplace_y 

        inertio = np.array(
            [
                xi*xi
                for xi in self.spl_x.collocPoints
                for yi in self.spl_y.collocPoints
            ]
        )
        self.inertio = sp.sparse.spdiags(inertio, 0, inertio.size, inertio.size)


    def init_potential(self, reg_x=0.03, reg_y=0.25):
        self.reg_x = reg_x
        self.reg_y = reg_y
        m1, m2, m3 = self.m1, self.m2, self.m3
        self.lv = lambda x, y: 1.0 / math.sqrt(x * x + self.reg_x) - (1.0/ math.sqrt((y - x / (1.0 + m2 / m1)) * (y - x / (1.0 + m2 / m1)) + self.reg_y)+ 1.0/ math.sqrt((y + x / (1.0 + m2 / m1)) * (y + x / (1.0 + m2 / m1)) + self.reg_y))

        potential = np.array(
            [
                self.lv(xi, yi)
                for xi in self.spl_x.collocPoints
                for yi in self.spl_y.collocPoints
            ]
        )
        self.v = sp.sparse.spdiags(potential, 0, potential.size, potential.size)
        self.v = 0.5 * self.v @ self.eye


    def init_hamiltonian(self):
        self.h = self.laplace + self.v


    def splinef(self, evec, x=1.0, y=0.0):
        values = np.array(
            [
               self.spl_x.fBSplineBC(x, i)*self.spl_y.fBSplineBC(y, j)
                for i in range(self.spl_x.getSpaceDim())
                for j in range(self.spl_y.getSpaceDim())
            ]
        )
        return evec.dot(values)


    def init_absorption(self, y0=1.0, k=1.0):
        absorb = lambda x, y: 0.0 if abs(y) < y0 else -1.0j * k * math.pow(y - y0, 2)
        
        absorption = np.array(
            [
                absorb(xi, yi)
                for xi in self.spl_x.collocPoints
                for yi in self.spl_y.collocPoints
            ]
        )
        self.absorb = sp.sparse.spdiags(absorption, 0, absorption.size, absorption.size)
        self.absorb = self.absorb @ self.eye



    def scale(self, t, t0, v):
        if t >= t0:
            Ra2 = 1.0 + math.pow(v*(t - t0), 4)
            Ra = math.sqrt(math.sqrt(Ra2))
            d2Ra = 3*math.pow(v, 4)*math.pow(t - t0, 2)/math.pow(Ra, 7)

            self.laplace = 1.0/(Ra*Ra)*self.laplace_x
            self.laplace += self.laplace_y
            self.laplace += 0.5*self.mu_y*Ra*d2Ra*self.inertio@self.eye


            potential = np.array(
                [
                    self.lv(Ra*xi, yi)
                    for xi in self.spl_x.collocPoints
                    for yi in self.spl_y.collocPoints
                ]
            )

            self.v = sp.sparse.spdiags(potential, 0, potential.size, potential.size)
            self.v = 0.5*self.v@self.eye
            

    def exp_hamiltonian(self, dt=1e-3, **kwargs):
        forward = self.eye  - 0.5 * dt * self.h
        backward = self.eye + 0.5 * dt * self.h
        
        # mv = lambda x: backward @ x
        # op = LinearOperator(backward.shape, matvec=mv)
        print('init')
        solve = factorized(backward)
        mv = lambda x: solve(forward @ x)
        op = LinearOperator(backward.shape, matvec=mv)

        evals, evecs = eigs(op, **kwargs)
        evals = -np.log(evals) / dt
        for i in range(len(evals)):
            norm=sp.dot(evecs[:,i], self.nMatr @ evecs[:,i])
            evecs[:,i]=evecs[:,i]/math.sqrt(norm)
            phasenorm=self.splinef(evecs[:,i], 1.0, 0.0)
            phasenorm=phasenorm/abs(phasenorm)
            evecs[:,i]=evecs[:,i]/phasenorm
        self.eigenvalues, self.eigenvectors=evals, evecs
        return evals, evecs

    def evolution_step(self, evec, dt=1e-3, iters=100, debug=False):
        v_forward = (self.eye  - 0.25j * dt * self.v)
        v_backward = (self.eye + 0.25j * dt * self.v)
        v_factorized = factorized(v_backward)
        pot_mv = lambda x: v_factorized(v_forward @ x)
        pot_op = LinearOperator(v_backward.shape, matvec=pot_mv)


        lpl_forward = self.eye  - 0.5j * dt * (self.laplace)
        lpl_backward = self.eye + 0.5j * dt * (self.laplace)
        print(lpl_forward.todense())
        lpl_factorized = factorized(lpl_backward)
        lpl_mv = lambda x: lpl_factorized(lpl_forward @ x)


        h_forward = (self.eye  - 0.5j * dt * (self.h))
        h_backward = (self.eye + 0.5j * dt * (self.h))
        print(h_forward.todense())
        h_factorized = factorized(h_backward)
        h_mv = lambda x: h_factorized(h_forward @ x)

        cache = []
        for i in range(iters):
            if debug:
                print(np.linalg.norm(evec.conjugate().dot(self.nMatr @ evec)))  
            evec = h_mv(evec)          
            #evec = pot_op(evec)
            #evec = lpl_mv(evec)
            #evec = pot_op(evec)
            cache.append(evec)

        #print(evec.dot(self.nMatr @ evec))
        return evec, cache


    def old_evol(self, evec, dt=1e-3,t0_scale=0.0, iters=100, debug=False):
        cache = []
        t = 0
        for i in range(iters):
            if t >= t0_scale:
                self.scale(t+dt*0.5, t0_scale, 0.0025)
            forward = self.eye  - 0.5j * dt * (self.h + self.absorb)
            backward = self.eye + 0.5j * dt * (self.h + self.absorb) 

            mv = lambda x: spsolve(backward, forward @ x)
            op = LinearOperator(backward.shape, matvec=mv)


            if debug:
                print(np.linalg.norm(evec.conjugate().dot(self.nMatr @ evec)))
            evec = op(evec)

            cache.append(evec.copy())

        return evec, cache
