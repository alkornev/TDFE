#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif
#include "hamiltonian3d.h"
#include "s32.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/Eigenvalues> 
#include <Eigen/QR>
#include <unsupported/Eigen/KroneckerProduct>

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <cmath>
#include <array>
#include <iostream>
#include <complex>
#include <typeinfo>


//double absorp_potential(double axi, double bxi):
    

Hamiltonian3D::Hamiltonian3D(
    const Eigen::Ref<const RVector>& aGrid, 
    const Eigen::Ref<const RVector>& bGrid, 
    const std::array<int, 4>& BCs,
    const std::array<double, 3>& masses, 
    double regularization, 
    int n
    ) 
: aSplines(aGrid, BCs[0], BCs[1]), bSplines(bGrid, BCs[2], BCs[3])
{
    r = regularization;
    // the last value in the array is the mass of electron
    m = masses;//{1836.0, 3671.0, 1.0};
    signs = {1.0, -1.0, -1.0};
    for (int k=0; k < (n % 3); k++){
        std::next_permutation(signs.begin(), signs.end());
        std::next_permutation(m.begin(), m.end());
        //std::cout << signs[0] << " " << signs[1] << " " << signs[2] << std::endl;
    }
    

    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    I = Eigen::SparseMatrix<double>(aNMax, aNMax);
    J = Eigen::SparseMatrix<double>(bNMax, bNMax);
    I.setIdentity();
    J.setIdentity();

    std::cout << "Initialization..." << std::endl;
    nMatr = generateNMatr();

    std::cout << "Normalization Matrix is ready! Continuing..." << std::endl;
    pMatr = generatePMatr();

    std::cout << "Collocation Matrix is ready! Continuing..." << std::endl;
    
    pMatrInv = generateInvPMatr();
    std::cout << "Inversion of Collocation Matrix is ready! Continuing..." << std::endl;

    h = generateTheHamiltonian();
    std::cout << "Hamiltonian is ready! Continuing..." << std::endl;
    //std::cout << h << std::endl;
}

Hamiltonian3D::~Hamiltonian3D(){}

SparseRMatrix Hamiltonian3D::generateNMatr()
{
    anMatr = aSplines.generateSNMatr();
    bnMatr = bSplines.generateSNMatr();
    SparseRMatrix nMatrix = Eigen::kroneckerProduct(anMatr, bnMatr).eval();

    return nMatrix; 
}

SparseRMatrix Hamiltonian3D::generatePMatr()
{
    apMatr = aSplines.generateSPMatr();
    bpMatr = bSplines.generateSPMatr();

    
    SparseRMatrix pMatrix = Eigen::kroneckerProduct(apMatr, bpMatr).eval();
    
    return pMatrix;   
}

SparseRMatrix Hamiltonian3D::generateInvPMatr()
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(apMatr);
    apMatrInv = solver.transpose().solve(I).eval();
    apMatrInv = apMatrInv.transpose();

    solver.compute(bpMatr);
    bpMatrInv = solver.transpose().solve(J).eval();
    bpMatrInv = bpMatrInv.transpose();
    

    SparseRMatrix pMatrixInv = Eigen::kroneckerProduct(apMatrInv, bpMatrInv).eval();
    
    return pMatrixInv;   
}

SparseRMatrix Hamiltonian3D::generateTheHamiltonian()
{
    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;
    //std::cout << "size of matrix: " << nMax << std::endl;
    SparseRMatrix ham(nMax, nMax);
    ham.reserve(Eigen::VectorXi::Constant(nMax,nMax));
    double potential = 0.0;

    double mu = 2.0*m[0]*m[1] / (m[0] + m[1]);
    double mu123 = 2.0*((m[0]+m[1])*m[2])/(m[0] + m[1] + m[2]);

    aLaplace = SparseRMatrix(aNMax, aNMax);
    bLaplace = SparseRMatrix(bNMax, bNMax);

    vxy = SparseRMatrix(nMax, nMax);

    for(int i = 0; i < aNMax; i++){
        double axi = aSplines.space.collocGrid[i];
        for(int j = 0; j < aNMax; j++){
            if (aSplines.d2BSplineBC(axi, j)) {
                aLaplace.insert(i,j) = -1.0/mu*aSplines.d2BSplineBC(axi, j);
            }
        }
    }

    for(int i = 0; i < bNMax; i++){
        double bxi = bSplines.space.collocGrid[i];
        for(int j = 0; j < bNMax; j++){
            if (bSplines.d2BSplineBC(bxi, j) != 0.0) {
                bLaplace.insert(i,j) = -1.0/mu123*bSplines.d2BSplineBC(bxi, j);
            }
        }
    }
    aLaplace.makeCompressed();
    bLaplace.makeCompressed();
    
    SparseRMatrix laplacian = Eigen::kroneckerProduct(aLaplace, bpMatr) + Eigen::kroneckerProduct(apMatr, bLaplace);
    laplacian.makeCompressed();
    std::cout << "laplacia non-zeros " << laplacian.nonZeros() << std::endl;

    double basisElement = 0.0;
    for (int i = 0; i < nMax; i++){
        //std::cout << "size of matrix: " << i << " " << std::endl;
        double axi = aSplines.space.collocGrid[i / bNMax];
        double bxi = bSplines.space.collocGrid[i % bNMax];
        for (int j = 0; j < nMax; j++){
            basisElement = aSplines.fBSplineBC(axi, j / bNMax)*bSplines.fBSplineBC(bxi, j % bNMax);
            if (basisElement != 0.0){
                potential = signs[0]/sqrt(axi*axi + 0.03) +
                        signs[1]/sqrt((bxi - m[0]*axi/(m[0]+m[1]))*(bxi - m[0]*axi/(m[0]+m[1])) + r) +
                        signs[2]/sqrt((bxi + m[1]*axi/(m[0]+m[1]))*(bxi + m[1]*axi/(m[0]+m[1])) + r);
                ham.insert(i, j) = 0.5*potential*basisElement;
            }
        }
    }
    ham = ham + laplacian;
    //std::cout << "Hamiltonian almost ready!" << std::endl;
    //std::cout << pMatrInv << std::endl;
    //std::cout << ham << std::endl;
    ham.makeCompressed();
    //std::cout << ham << std::endl;

    return ham;   
}

double Hamiltonian3D::getEigenfunction(
    const Eigen::Ref<const RVector>& coefs,
    double x,
    double y
)
{
    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;

    //double t = 0.0;
    Eigen::VectorXd t = Eigen::VectorXd(nMax);
    //std::cout <<"size of points: "<< t.size() << std::endl;

    //int aStartIdx = aSplines.locate(x);
    //int bStartIdx = bSplines.locate(y);

    int idx = 0;
    
    for(int i = 0; i < nMax; i++){
        t[i] = aSplines.fBSplineBC(x, i / bNMax)*bSplines.fBSplineBC(y, i % bNMax);
    }
    // int aShift = (aSplines.rightBC >= 0 && aStartIdx == (aSplines.splineBCdim - 1)) || (aSplines.leftBC >= 0 && aStartIdx == 0);
    // int bShift = (bSplines.rightBC >= 0 && bStartIdx == (bSplines.splineBCdim - 1)) || (bSplines.leftBC >= 0 && bStartIdx == 0);

    // for(int i = aStartIdx; i < aStartIdx + aSplines.space.splinesPerNode * 2; i++){
    //     for (int j = bStartIdx; j < bStartIdx + bSplines.space.splinesPerNode * 2; j++){
    //         idx = j + bNMax*i;
    //         t += coefs(idx) * aSplines.fBSplineBC(x, idx / bNMax)*bSplines.fBSplineBC(y, idx % bNMax);
    //     }
    // }
    double res = t.dot(coefs);

    return res;
}

void Hamiltonian3D::getTheSpectrum(int vector_n, int krylov_n)
{   
    auto hsinv = h;
    auto test = Eigen::kroneckerProduct(I, bpMatrInv);

    std::cout << "h" << h.nonZeros() << std::endl;
    std::cout << "apMatr" << apMatr.nonZeros() << std::endl;
    std::cout << "bpMatr" << bpMatr.nonZeros() << std::endl;
    std::cout << "aLaplace" << aLaplace.nonZeros() << std::endl;
    std::cout << "bLaplace" << bLaplace.nonZeros() << std::endl;


    
    hsinv = Eigen::kroneckerProduct(I, bpMatrInv)*hsinv;
    hsinv = Eigen::kroneckerProduct(apMatrInv, J)*hsinv;
    std::cout << "hsinv" << hsinv.nonZeros() << std::endl;

    Spectra::SparseGenMatProd<double, Eigen::RowMajor> op(hsinv);

    // Construct eigen solver object, requesting the largest three eigenvalues
    Spectra::GenEigsSolver<Spectra::SparseGenMatProd<double, Eigen::RowMajor>> eigs(op, vector_n, krylov_n);
    eigs.init();

    int nconv = eigs.compute(Spectra::SortRule::SmallestReal);
    
    this->evectors = eigs.eigenvectors().real();
    this->evalues = eigs.eigenvalues();

    Eigen::VectorXcd evalues;
    if (eigs.info() == Spectra::CompInfo::Successful){
        evalues = eigs.eigenvalues();
        std::cout << evalues << std::endl;
    } else {
        std::cout << "Error" << std::endl;
    }
        

    for(int i=0; i<evalues.size(); i++){
        double norm = sqrt(evectors.col(i).dot(nMatr*evectors.col(i)));
        evectors.col(i) /=norm;
        double phasenorm = getEigenfunction(this->evectors.col(i), 1.0, 1.0);
        phasenorm = phasenorm/abs(phasenorm);
        evectors.col(i) /= phasenorm;
    }
}

CVector Hamiltonian3D::getEigenvalues(){ return evalues; }

RMatrix Hamiltonian3D::getEigenvectors(){ return evectors; }

SparseRMatrix Hamiltonian3D::getPMatr(){ return pMatr; }

SparseRMatrix Hamiltonian3D::getNMatr(){ return nMatr; }

SparseRMatrix Hamiltonian3D::getHamiltonian(){ return h; }

SparseCMatrix Hamiltonian3D::getImpulse(double t)
{
    double V0 = std::sqrt(I0/3.5095);
    double envArg = (t + t0)/tau;

    int bNMax = bSplines.splineBCdim;
    SparseCMatrix eField = SparseCMatrix(bNMax, bNMax);
    
    eField.reserve(Eigen::VectorXi::Constant(bNMax,bNMax));
    double bxi = 0.0;

    double pulse = V0*std::exp(-envArg*envArg)*cos(w0*(t + t0) + phi);
    
        std::cout << V0 << " " << pulse << " " << w0 << " " << phi << std::endl;

    for(int i=0; i < bNMax; i++){
        bxi = bSplines.space.collocGrid[i];
        for (int j=0; j < bNMax; j++){
            if (bpMatr.coeff(i, j) != 0.0) {
                eField.insert(i, j) = -(1.0+m[2]/(m[0]+m[1]+m[2]))*bxi*pulse*bpMatr.coeff(i, j);
            }
        }
    }
    eField.makeCompressed();

    // std::cout << "eField after filling\n" << eField << std::endl;
    // std::cout << "bpMatrInv\n" << bpMatrInv << std::endl;
    // std::cout << "bpMatr\n" << bpMatr << std::endl;

    //eField = bpMatrInv * eField;

    //electricField = Eigen::kroneckerProduct(apMatr, eField);
    // electricField = Eigen::kroneckerProduct(I, bpMatrInv)*electricField;
    // electricField = Eigen::kroneckerProduct(apMatrInv, J)*electricField;
    //electricField = Eigen::kroneckerProduct(I, eField);
    //std::cout << eField << std::endl;
    //
    //std::cout << electricField << std::endl;
    //electricField.makeCompressed();
    return eField;
}

CVector Hamiltonian3D::evolutionStep(CVector state, int iter)
{
    /*
        We want to calculate this formula:
        inv(I + i/2*A)*(I - i/2*A) x = c
        decompose it into two steps:
        first step:
        (I - i/2*A) y = c
        second step:
        inv(I + i/2*A) x = y 
        second step is equal to the following system:
        (I + i/2*A) y = x
    */
    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;

    std::complex<double> im(0.0, dt);
    CVector y;
    
    y = S*state;
    y = solver.solve(y);

    //std::cout << "condition " << t0 + dt*iter << std::endl;

    Eigen::Map<Eigen::MatrixXcd> map(y.data(), bNMax, aNMax); 
    std::cout << "t " << t0 + dt*iter << ", iter " << iter;

    if (std::abs(t0 + dt*iter) < tau) {

        //std::cout << dt*iter << std::endl;
        auto Vt = getImpulse(dt*(iter + 0.5));
        //std::cout << Vt << std::endl;
        Eigen::SparseMatrix<std::complex<double>> St = bpMatr.cast<std::complex<double>>() - 0.5 * im * (Vt + W);
        Eigen::SparseMatrix<std::complex<double>> Ft = bpMatr.cast<std::complex<double>>() + 0.5 * im * (Vt + W);

        Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solverVt;
        // Eigen::SparseMatrix<std::complex<double>> FKt = Eigen::kroneckerProduct(myJ, Ft);
        // Eigen::SparseMatrix<std::complex<double>> SKt = Eigen::kroneckerProduct(myJ, St);
        solverVt.compute(Ft);

        map = St*map;
        map = solverVt.solve(map);
    }
    
    y = S*y;    
    y = solver.solve(y);

    std::cout << ", norm " << std::sqrt(y.dot(nMatr*y)) << std::endl;

    return y;
}

void Hamiltonian3D::initHamiltonianLU(){
    std::complex<double> im(0.0, dt);

    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;
    eye = SparseCMatrix(nMax, nMax);
    eye.setIdentity();

    S = pMatr.cast<std::complex<double>>() - 0.25 * im * h.cast<std::complex<double>>();
    Eigen::SparseMatrix<std::complex<double>> F = pMatr.cast<std::complex<double>>() + 0.25 * im * h.cast<std::complex<double>>();

    solver.compute(F);
}

void Hamiltonian3D::initImpulse(double init_time, double phase, double intensity, double freq, double duration, double step){
    w0 = freq;
    I0 = intensity;
    phi = phase;
    tau = duration;
    t0 = init_time;
    dt = step;
}

void Hamiltonian3D::initAbsorption(double alpha, double width){
    std::complex<double> im(0.0, 1.0/dt);
    int bNMax = bSplines.splineBCdim;
    W = SparseCMatrix(bNMax, bNMax);
    
    W.reserve(Eigen::VectorXi::Constant(bNMax,bNMax));
    double bxi = 0.0;
    double R = bSplines.space.grid[bSplines.space.grid.size()-1];

    //std::cout << R << std::endl;
    for(int i=0; i < bNMax; i++){
        bxi = bSplines.space.collocGrid[i];
        for (int j=0; j < bNMax; j++){
            //std::cout << "bpmatr" << bpMatr.coeff(i, j) << " bxi " << bxi << "pulse " << pulse;
            if (bpMatr.coeff(i, j) != 0.0) {
                W.insert(i, j) = im*std::log(1.0 - std::pow(std::cos(M_PI/2.0*(1.0 - std::abs(bxi)/R)), 50))*bpMatr.coeff(i, j);
            }
        }
    }
    W.makeCompressed();
    //std::cout << W << std::endl;
}
