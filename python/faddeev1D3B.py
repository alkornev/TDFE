from distutils.fancy_getopt import fancy_getopt
import math
import numpy as np
import scipy as sp
import scipy.integrate
import scipy.sparse
import scipy.io
from itertools import product
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
from python.hamiltonian1D3B import Hamiltonian1D3B


class Component:
    def __init__(self, masses, regulizers, pair_index, grid_x, grid_y, left_bcx=1, right_bcx=1, left_bcy=1, right_bcy=1):
        self.masses = masses
        self.regulizers = regulizers
        self.signs = [-1, -1, 1]
        self.grid_x = grid_x
        self.grid_y = grid_y
        self.format = 'csc'
        self.spl_x = PyHermiteBC(grid_x, left_bcx, right_bcx)
        self.spl_y = PyHermiteBC(grid_y, left_bcy, right_bcy)

        self.ident_x = self.spl_x.getSPMatr()
        self.ident_y = self.spl_y.getSPMatr()
        
        self.Ix = sp.sparse.eye(self.spl_x.getSpaceDim(), format=self.format)
        self.Iy = sp.sparse.eye(self.spl_y.getSpaceDim(), format=self.format)

        self.idx = pair_index
        self.nx = self.spl_x.getSNMatr()
        self.ny = self.spl_y.getSNMatr()

        self.nMatr = sp.sparse.kron(self.nx, self.ny, format=self.format)
        self.ident = sp.sparse.kron(self.ident_x, self.ident_y, format=self.format)

        self.init_laplacian()
        self.init_potential()
        self.dim = self.spl_x.getSpaceDim() * self.spl_y.getSpaceDim()
        self.decompose()


    def init_laplacian(self):
        alpha = self.idx % 3
        beta = (self.idx + 1) % 3
        gamma = (self.idx + 2) % 3
        self.mu_1 = 0.5 * (1/self.masses[beta] + 1/self.masses[gamma])
        self.mu_2 = 0.5 * (1/(self.masses[beta] + self.masses[gamma]) + 1/self.masses[alpha])
        
        self.ddx = np.array(
            [
                [self.spl_x.d2BSplineBC(xi, i) for i in range(self.spl_x.getSpaceDim())]
                for j, xi in enumerate(self.spl_x.collocPoints)
            ]
        )
        self.ddy = np.array(
            [
                [self.spl_y.d2BSplineBC(yi, i) for i in range(self.spl_y.getSpaceDim())]
                for j, yi in enumerate(self.spl_y.collocPoints)
            ]
        )

        self.ddx = sp.sparse.csc_matrix(self.ddx)
        self.ddy = sp.sparse.csc_matrix(self.ddy)

        self.laplace_x = -self.mu_1 * sp.sparse.kron(self.ddx, self.ident_y, format=self.format)
        self.laplace_y = -self.mu_2 * sp.sparse.kron(self.ident_x, self.ddy, format=self.format)
        self.laplace = self.laplace_x + self.laplace_y

    def decompose(self):
        Dx = -self.mu_1 * self.ddx.toarray()
        Sx = self.ident_x.toarray()
        Ax = (Dx + self.vx @ Sx)

        Dy = -self.mu_2 * self.ddy.toarray()
        Sy = self.ident_y.toarray()
        print(f"reduced mass: m1={1/self.mu_1}, m2={1/self.mu_2}")


        xvals, Wxl, Wxr = sp.linalg.eig(Ax, b=Sx,  left=True, right=True)
        Sx_diag = np.diag(Wxl.T @ Sx @ Wxr)
        print(f"index {self.idx}")
        Sx_sqrt_diag_inv = sp.sparse.diags(1/np.sqrt(Sx_diag), format=self.format)

        self.Wxl = sp.sparse.csc_matrix(Sx_sqrt_diag_inv @ Wxl.T)
        self.Wxr = sp.sparse.csc_matrix(Wxr @ Sx_sqrt_diag_inv)
        xvals = sp.sparse.diags(xvals, format=self.format)
        self.xvals = sp.sparse.csc_matrix(xvals)
        yvals, Wyl, Wyr = sp.linalg.eig(Sy, b=Dy,  left=True, right=True)
        Dy_diag = np.diag(Wyl.T @ Dy @ Wyr)
        print("yvals")
        Dy_sqrt_diag_inv = sp.sparse.diags(1/np.sqrt(Dy_diag), format=self.format)

        self.Wyl = sp.sparse.csc_matrix(Dy_sqrt_diag_inv @ Wyl.T)
        self.Wyr = sp.sparse.csc_matrix(Wyr @ Dy_sqrt_diag_inv)
        yvals = sp.sparse.spdiags(yvals, 0, yvals.size, yvals.size, format=self.format)

        self.yvals = sp.sparse.csc_matrix(yvals)
        self.Wr = lambda x: sp.sparse.kron(self.Wxl, self.Iy) @ (sp.sparse.kron(self.Ix, self.Wyl) @ x)
        self.Wl = lambda x: sp.sparse.kron(self.Wxr, self.Iy) @ (sp.sparse.kron(self.Ix, self.Wyr) @ x)

    def diag_inverse(self, E):
        self.diag = (sp.sparse.kron(self.xvals, self.yvals) + sp.sparse.kron(self.Ix, self.Iy) - E * sp.sparse.kron(self.Ix, self.yvals)).toarray()
        self.diag_inv = sp.sparse.diags(1.0/np.diag(self.diag), format=self.format)
        return self.diag_inv


    def init_potential(self):
        reg = self.regulizers[self.idx]
        self.lv = lambda x, y: 0.5*self.signs[self.idx] / math.sqrt(x * x + reg)

        potential = np.array(
            [
                self.lv(xi, yi)
                for xi in self.spl_x.collocPoints
                for yi in self.spl_y.collocPoints
            ]
        )
        vx = np.array(
            [
                self.lv(xi, 0) for xi in self.spl_x.collocPoints
            ]
        )
        self.vx = sp.sparse.spdiags(vx, 0, vx.size, vx.size, format=self.format)
        self.v = sp.sparse.spdiags(potential, 0, potential.size, potential.size, format=self.format)
        self.v =  self.v

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


class MulticomponentRHS:
    def __init__(self, masses):
        self.transforms = np.zeros((3, 3, 2, 2))
        self.masses = masses
        total_mass = sum(masses)
        print(masses)
        print(total_mass)
        # i - beta
        # j - alpha
        dets = []
        for alpha in range(3):
            for beta in range(3):
                if alpha != beta:
                    pwr = (-1)**(-alpha+beta)
                    sgn = (alpha-beta)/abs(alpha-beta)
                    c1 = -masses[alpha]/(total_mass - masses[beta])
                    c2 = -masses[beta]/(total_mass - masses[alpha])
                    s1 = -sgn*pwr
                    s2 = pwr*sgn*((total_mass - masses[alpha] - masses[beta])*total_mass/(total_mass - masses[alpha])/(total_mass - masses[beta]))
                else:
                    c1 = c2 = 1.0
                    s1 = s2 = 0.0
                
                    # pwr = (-1)**(-alpha+beta)
                    # sgn = (alpha-beta)/abs(alpha-beta)
                    # c = -math.sqrt(masses[alpha]*masses[beta]/(M - masses[alpha])/(M - masses[beta]))
                    # s = -pwr*sgn*math.sqrt(M*(M - masses[alpha] -masses[beta])/(M - masses[alpha])/(M - masses[beta]))
                    # lx = 1/(0.5 * (1/masses[beta] + 1/(M-masses[alpha]-masses[beta])))
                    # ly = 1/(0.5 * (1/(M-masses[alpha]) + 1/masses[alpha]))
                    # rx = 0.5 * (1/masses[alpha] + 1/(M-masses[alpha]-masses[beta]))
                    # ry = 0.5 * (1/(M-masses[beta]) + 1/masses[beta])


                self.transforms[beta][alpha][0][0] = c1
                self.transforms[beta][alpha][1][1] = c2
                self.transforms[beta][alpha][0][1] = s1
                self.transforms[beta][alpha][1][0] = s2

                #print(f"components {i}, {j}, {np.linalg.det(self.transforms[i][j])}")
                dets.append(np.linalg.det(self.transforms[beta][alpha]))
        print(dets)

    def convert_gridpoints(self, to: Component, frm: Component):
        orig_colloc_grids = np.array([pair for pair in product(frm.spl_x.collocPoints, frm.spl_y.collocPoints)])
        new_colloc_grids = self.transforms[frm.idx][to.idx] @ orig_colloc_grids.T
        return new_colloc_grids

    def transform(self, to: Component, frm: Component):
        orig_colloc_grids = np.array([pair for pair in product(frm.spl_x.collocPoints, frm.spl_y.collocPoints)])
        new_colloc_grids = self.transforms[frm.idx][to.idx] @ orig_colloc_grids.T

        rhs_row_dim = frm.spl_x.getSpaceDim()*frm.spl_y.getSpaceDim()
        rhs_col_dim = to.spl_x.getSpaceDim()*to.spl_y.getSpaceDim()
        rhs = np.zeros((rhs_row_dim, rhs_col_dim))

        for idx, (xi, yj) in enumerate(zip(*new_colloc_grids)):
            row_i = idx // to.spl_y.getSpaceDim()
            row_j = idx % to.spl_y.getSpaceDim()
            row_idx = row_i*to.spl_y.getSpaceDim() + row_j
            for i in range(to.spl_x.getSpaceDim()):
                if to.spl_x.fBSplineBC(xi, i) == 0.0:
                    continue
                for j in range(to.spl_y.getSpaceDim()):
                    if to.spl_y.fBSplineBC(yj, j) == 0.0:
                        continue
                    col_idx = i*to.spl_y.getSpaceDim() + j
                    rhs[row_idx, col_idx] = to.spl_x.fBSplineBC(xi, i)*to.spl_y.fBSplineBC(yj, j)

        rhs_ident = sp.sparse.csc_matrix(rhs) 
        return rhs_ident



class Faddeev1D3B:
    def __init__(self, grids, masses=None,regulizers=None, left_bc=1, right_bc=1):
        if masses is None:
            self.m1 = 1836.0
            self.m2 = 1836.0
            self.m3 = 1.0
            self.masses = [self.m1, self.m2, self.m3]
        else:
            self.m1 = masses[0]
            self.m2 = masses[1]
            self.m3 = masses[2]
            self.masses = [self.m1, self.m2, self.m3]

        self.format = 'csc'
        if regulizers is None:
            self.regulizers = [0.25, 0.25, 0.03]
        else:
            self.regulizers = regulizers

        self.alpha = Component(self.masses, self.regulizers, 0, grids[0], grids[1])
        self.beta = Component(self.masses, self.regulizers, 1, grids[2], grids[3])
        self.gamma = Component(self.masses, self.regulizers, 2, grids[4] ,grids[5])
        
        self.nMatr = sp.sparse.bmat(([[self.alpha.nMatr, None, None], [None, self.beta.nMatr, None], [None, None, self.gamma.nMatr]]), format=self.format)
    
    def preconditioner(self, E):
        def mv(t):
            t_alpha = t[:self.alpha.dim]
            t[:self.alpha.dim] = self.alpha.Wl(self.alpha.diag_inverse(E) @ self.alpha.Wr(t_alpha))
            beta_bound = self.alpha.dim+self.beta.dim
            t_beta = t[self.alpha.dim:beta_bound]
            t[self.alpha.dim:beta_bound] = self.beta.Wl(self.beta.diag_inverse(E) @ self.beta.Wr(t_beta))
            t_gamma = t[beta_bound:]
            t[beta_bound:] = self.gamma.Wl(self.gamma.diag_inverse(E) @ self.gamma.Wr(t_gamma))
            t = self.c @ t
            t = -self.v @ t
            return t
        
        return mv

    def init_fhamiltonian(self):
        self.h0 = sp.sparse.bmat(([[self.alpha.laplace, None, None], [None, self.beta.laplace, None], [None, None, self.gamma.laplace]]), format=self.format)
        self.ident = sp.sparse.bmat(([[self.alpha.ident, None, None], [None, self.beta.ident, None], [None, None, self.gamma.ident]]), format=self.format)
        self.v = sp.sparse.bmat(([[self.alpha.v, None, None], [None, self.beta.v, None], [None, None, self.gamma.v]]), format=self.format)
        rhs = self.rhs = MulticomponentRHS(self.masses)
        from_alpha = self.from_alpha = [None, rhs.transform(self.beta, self.alpha), rhs.transform(self.gamma, self.alpha)]
        from_beta = self.from_beta = [rhs.transform(self.alpha, self.beta), None, rhs.transform(self.gamma, self.beta)]
        from_gamma = self.from_gamma = [rhs.transform(self.alpha, self.gamma), rhs.transform(self.beta, self.gamma), None]

        self.c = sp.sparse.bmat([from_alpha, from_beta, from_gamma], format=self.format)
        self.vc = self.v @ self.c

        self.h = self.h0 + self.v @ self.ident + self.v @ self.c

    def splinef(self, evec, x, y, idx=0):
        alpha_dim = self.alpha.dim
        beta_dim = self.beta.dim
        gamma_dim = self.gamma.dim
        values = np.zeros(alpha_dim + beta_dim + gamma_dim)

        alpha_xy = self.rhs.transforms[idx][self.alpha.idx] @ np.array([[x], [y]])
        xa,ya = alpha_xy[0], alpha_xy[1]

        for i in range(self.alpha.spl_x.getSpaceDim()):
            if self.alpha.spl_x.fBSplineBC(xa, i) == 0.0:
               continue
            for j in range(self.alpha.spl_y.getSpaceDim()):
                if self.alpha.spl_y.fBSplineBC(ya, j) == 0.0:
                   continue
                values[i*self.alpha.spl_y.getSpaceDim() + j] = self.alpha.spl_x.fBSplineBC(xa, i) * self.alpha.spl_y.fBSplineBC(ya, j)

        beta_xy = self.rhs.transforms[idx][self.beta.idx] @ np.array([[x], [y]])
        xb,yb = beta_xy[0], beta_xy[1]
        for i in range(self.beta.spl_x.getSpaceDim()):
            if self.beta.spl_x.fBSplineBC(xb, i) == 0.0:
               continue
            for j in range(self.beta.spl_y.getSpaceDim()):
                if self.beta.spl_y.fBSplineBC(yb, j) == 0.0:
                   continue
                values[alpha_dim + i*self.beta.spl_y.getSpaceDim() + j] = self.beta.spl_x.fBSplineBC(xb, i) * self.beta.spl_y.fBSplineBC(yb, j)

        gamma_xy = self.rhs.transforms[idx][self.gamma.idx] @ np.array([[x], [y]])
        xg,yg = gamma_xy[0], gamma_xy[1]

        for i in range(self.gamma.spl_x.getSpaceDim()):
            if self.gamma.spl_x.fBSplineBC(xg, i) == 0.0:
               continue
            for j in range(self.gamma.spl_y.getSpaceDim()):
                if self.gamma.spl_y.fBSplineBC(yg, j) == 0.0:
                   continue
                values[alpha_dim + beta_dim + i*self.gamma.spl_y.getSpaceDim() + j] = self.gamma.spl_x.fBSplineBC(xg, i) * self.gamma.spl_y.fBSplineBC(yg, j)


        return evec.dot(values)

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
                norm = sp.dot(evecs[:, i], self.nMatr @ evecs[:, i])
                evecs[:, i] = evecs[:, i] / math.sqrt(norm)
                phasenorm = self.splinef(evecs[:, i], 1.0, 0.0, idx=2)
                phasenorm = phasenorm / abs(phasenorm)
                evecs[:, i] = evecs[:, i] / phasenorm
            self.eigenvalues, self.eigenvectors = evals, evecs
            return evals, evecs

    def precondition_hamiltonian(self, lE=-0.7, rE=-0.6, debug=False, **kwargs):
            #forward = self.ident - 0.5 * dt * self.h
            #backward = self.ident + 0.5 * dt * self.h

            lmbda = 1000
            while abs(lmbda - 1) > 1e-6: 
                E = 0.5 * (lE + rE)
                print(f"E: {E}")
                lop = LinearOperator(self.vc.shape, matvec=self.preconditioner(lE))
                lambda_l, _ = eigs(lop, **kwargs)
                print(f"left eigenvalue: {lambda_l}")
                rop = LinearOperator(self.vc.shape, matvec=self.preconditioner(rE))
                lambda_r, _ = eigs(rop, **kwargs)
                print(f"right eigenvalue: {lambda_r}")
                mop = LinearOperator(self.vc.shape, matvec=self.preconditioner(E))
                lmbda, evec =  eigs(mop, **kwargs)
                #print(f"E:   {E}, lmbda: {lmbda}, {lambda_l}, {lambda_r}")
                if np.sign(np.abs(lmbda) - 1) == np.sign(np.abs(lmbda) - 1):
                    rE = E
                else:
                    lE = E

            self.eigenvalues, self.eigenvectors = E, evec
            return E, evec


class ReducedFaddeev1D3B:
    def __init__(self, grids, masses=None,regulizers=None, left_bc=1, right_bc=1):
        if masses is None:
            self.m1 = 1836.0
            self.m2 = 1836.0
            self.m3 = 1.0
            self.masses = [self.m1, self.m2, self.m3]
        else:
            self.m1 = masses[0]
            self.m2 = masses[1]
            self.m3 = masses[2]
            self.masses = [self.m1, self.m2, self.m3]

        self.format = 'csc'
        if regulizers is None:
            self.regulizers = [0.25, 0.25, 0.03]
        else:
            self.regulizers = regulizers

        self.alpha = Component(self.masses, self.regulizers, 0, grids[0], grids[1])
        self.beta = Component(self.masses, self.regulizers, 1, grids[2], grids[3])
        self.gamma = Component(self.masses, self.regulizers, 2, grids[4], grids[5])

        self.nMatr = sp.sparse.bmat([[self.alpha.nMatr, None], [None, self.beta.nMatr]], format=self.format)

    def init_fhamiltonian(self):
        self.h0 = sp.sparse.bmat([[self.alpha.laplace, None], [None, self.beta.laplace]], format=self.format)
        self.ident = sp.sparse.bmat([[self.alpha.ident, None], [None, self.beta.ident]], format=self.format)
        self.v = sp.sparse.bmat([[self.alpha.v, None], [None, self.beta.v]], format=self.format)
        rhs = self.rhs = MulticomponentRHS(self.masses)
        from_alpha = self.from_alpha = [None, rhs.transform(self.beta, self.alpha)]
        from_beta = self.from_beta = [rhs.transform(self.alpha, self.beta), None]
        self.c = sp.sparse.bmat([from_alpha, from_beta], format=self.format)

        print(rhs.convert_gridpoints(self.gamma, self.alpha))
        v3_ga = [0.5/math.sqrt(x*x + 0.03) for x, y in zip(*rhs.convert_gridpoints(self.gamma, self.alpha))]
        v3_ga = sp.sparse.diags(v3_ga, format=self.format)

        v3_gb = [0.5/math.sqrt(x*x + 0.03) for x, y in zip(*rhs.convert_gridpoints(self.gamma, self.beta))]
        v3_gb = sp.sparse.diags(v3_gb, format=self.format)

        self.rhs_v = sp.sparse.bmat([[None, v3_ga], [v3_gb, None]], format=self.format)
        self.vc = self.rhs_v @ self.c

        self.h = self.h0 + (self.v + self.rhs_v)@ self.ident + self.v @ self.c 

    def splinef(self, evec, x, y, idx=0):
        alpha_dim = self.alpha.dim
        beta_dim = self.beta.dim
        values = np.zeros(alpha_dim + beta_dim)

        alpha_xy = self.rhs.transforms[idx][self.alpha.idx] @ np.array([[x], [y]])
        xa,ya = alpha_xy[0], alpha_xy[1]

        for i in range(self.alpha.spl_x.getSpaceDim()):
            if self.alpha.spl_x.fBSplineBC(xa, i) == 0.0:
               continue
            for j in range(self.alpha.spl_y.getSpaceDim()):
                if self.alpha.spl_y.fBSplineBC(ya, j) == 0.0:
                   continue
                values[i*self.alpha.spl_y.getSpaceDim() + j] = self.alpha.spl_x.fBSplineBC(xa, i) * self.alpha.spl_y.fBSplineBC(ya, j)

        beta_xy = self.rhs.transforms[idx][self.beta.idx] @ np.array([[x], [y]])
        xb,yb = beta_xy[0], beta_xy[1]
        for i in range(self.beta.spl_x.getSpaceDim()):
            if self.beta.spl_x.fBSplineBC(xb, i) == 0.0:
               continue
            for j in range(self.beta.spl_y.getSpaceDim()):
                if self.beta.spl_y.fBSplineBC(yb, j) == 0.0:
                   continue
                values[alpha_dim + i*self.beta.spl_y.getSpaceDim() + j] = self.beta.spl_x.fBSplineBC(xb, i) * self.beta.spl_y.fBSplineBC(yb, j)

        return evec.dot(values)

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
                norm = sp.dot(evecs[:, i], self.nMatr @ evecs[:, i])
                #evecs[:, i] = evecs[:, i] / math.sqrt(norm)
                #phasenorm = self.splinef(evecs[:, i], 1.0, 0.0, idx=2)
                #phasenorm = phasenorm / abs(phasenorm)
                #evecs[:, i] = evecs[:, i] / phasenorm
            self.eigenvalues, self.eigenvectors = evals, evecs
            return evals, evecs

    def preconditioner(self, E):
        def mv(t):
            t_alpha = t[:self.alpha.dim]
            t[:self.alpha.dim] = self.alpha.Wl(self.alpha.diag_inverse(E) @ self.alpha.Wr(t_alpha))
            beta_bound = self.alpha.dim+self.beta.dim
            t_beta = t[self.alpha.dim:beta_bound]
            t[self.alpha.dim:beta_bound] = self.beta.Wl(self.beta.diag_inverse(E) @ self.beta.Wr(t_beta))
            t_gamma = t[beta_bound:]
            t[beta_bound:] = self.gamma.Wl(self.gamma.diag_inverse(E) @ self.gamma.Wr(t_gamma))
            t = self.c @ t
            t = -self.v @ t
            return t
        
        return mv

    def precondition_hamiltonian(self, lE=-0.7, rE=-0.6, debug=False, **kwargs):
            #forward = self.ident - 0.5 * dt * self.h
            #backward = self.ident + 0.5 * dt * self.h
            lmbda = 1000
            while abs(lmbda - 1) > 1e-6: 
                E = 0.5 * (lE + rE)
                print(f"E: {E}")
                lop = LinearOperator(self.vc.shape, matvec=self.preconditioner(lE))
                lambda_l, _ = eigs(lop, **kwargs)
                print(f"left eigenvalue: {lambda_l}")
                rop = LinearOperator(self.vc.shape, matvec=self.preconditioner(rE))
                lambda_r, _ = eigs(rop, **kwargs)
                print(f"right eigenvalue: {lambda_r}")
                mop = LinearOperator(self.vc.shape, matvec=self.preconditioner(E))
                lmbda, evec =  eigs(mop, **kwargs)
                #print(f"E:   {E}, lmbda: {lmbda}, {lambda_l}, {lambda_r}")
                if np.sign(np.abs(lmbda) - 1) == np.sign(np.abs(lmbda) - 1):
                    rE = E
                else:
                    lE = E                    
            # for i in range(len(evals)):
            #     print(i)
            #     norm = sp.dot(evecs[:, i], self.nMatr @ evecs[:, i])
            #     evecs[:, i] = evecs[:, i] / math.sqrt(norm)
            #     phasenorm = self.splinef(evecs[:, i], 1.0, 0.0)
            #     phasenorm = phasenorm / abs(phasenorm)
            #     evecs[:, i] = evecs[:, i] / phasenorm
            self.eigenvalues, self.eigenvectors = E, evec
            return E, evec