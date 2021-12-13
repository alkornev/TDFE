#ifndef HAMILTONIAN3D_H
#define HAMILTONIAN3D_H

#include "s32.h"
#include "hamiltonian2d.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <array>
#include <complex>


typedef Eigen::VectorXd RVector;
typedef Eigen::VectorXcd CVector;
typedef Eigen::MatrixXd RMatrix;
typedef Eigen::MatrixXcd CMatrix;

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseRMatrix;
typedef Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> SparseCMatrix;


struct GroundState{
    double eigValue;
    RVector eigVector;
};


class Hamiltonian3D
{
public:
    double r;
    
    double w0;
    double I0;
    double phi;
    double tau;
    double t0;
    double dt;

    double scale_vel;
    double scale_t0;
    
    std::array<double, 3> m;  
    std::array<double, 3> signs;
    CHermiteBC aSplines;
    CHermiteBC bSplines;

    SparseRMatrix h;  
    SparseRMatrix aLaplace;
    SparseRMatrix bLaplace;
    SparseRMatrix vxy;

    SparseRMatrix apMatr;
    SparseRMatrix bpMatr;
    SparseRMatrix pMatr;

    SparseRMatrix apMatrInv;
    SparseRMatrix bpMatrInv;
    SparseRMatrix pMatrInv;

    SparseRMatrix anMatr;
    SparseRMatrix bnMatr;
    SparseRMatrix nMatr;

    RMatrix evectors;
    CVector evalues;

    SparseCMatrix eye;
    Eigen::SparseMatrix<double> I;
    Eigen::SparseMatrix<double> J;
    
    Eigen::SparseMatrix<std::complex<double>> S;
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;

    SparseCMatrix W;
    Hamiltonian3D(
        const Eigen::Ref<const RVector>& aGrid, 
        const Eigen::Ref<const RVector>& bGrid, 
        const std::array<int, 4>& BCs,
        const std::array<double, 3>& masses,
        double regularization, int n
        );
    ~Hamiltonian3D();
    SparseRMatrix generateNMatr();
    SparseRMatrix generatePMatr();
    SparseRMatrix generateInvPMatr();
    SparseRMatrix generateTheHamiltonian();
    double getEigenfunction(const Eigen::Ref<const RVector>& coefs, double x, double y);
    void getTheSpectrum(int vector_n, int krylov_n);
    void getTheSpectrumShiftInvert(int vector_n, int krylov_n);
    void initHamiltonianLU();
    void initScaling(double vel, double t0);
    void initAbsorption(int smooth);
    void initImpulse(double init_time, double init_phase, double intensity, double freq, double duration, double step);
    CVector getEigenvalues();
    RMatrix getEigenvectors();
    SparseRMatrix getPMatr();
    SparseRMatrix getNMatr();
    SparseRMatrix getHamiltonian();
    SparseCMatrix getImpulse(double t);
    void scalePotential(double t);
    void scaleLaplace(double t);
    GroundState getExpGroundState(double err, double dt);
    CVector evolutionStep(CVector state, int iter);
};


class ThreeBodyJacobiCoordinates
{
    std::array<double, 3> m;
public:
    ThreeBodyJacobiCoordinates(std::array<double, 3> m){
        // We arrange particles in the specific order.
        // ith == 1 means that 
    }
};
#endif