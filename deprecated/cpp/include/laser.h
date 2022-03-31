#ifndef LASER_H
#define LASER_H


#include "s32.h"
#include "hamiltonian2d.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <array>
#include <complex>


using RVector = Eigen::VectorXd;
using CVector = Eigen::VectorXcd;
using RMatrix = Eigen::MatrixXd;
using CMatrix = Eigen::MatrixXcd;
using SparseRMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using SparseCMatrix = Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>;


struct GroundState{
    double eigValue;
    RVector eigVector;
};


class LaserImpulse
{
public:
    double w0;
    double I0;
    double phi;
    double tau;
    double t0;
    double dt;

    CHermiteBC aSplines;
    CHermiteBC bSplines;

    SparseRMatrix apMatr;
    SparseRMatrix bpMatr;
    SparseRMatrix pMatr;

    SparseCMatrix W;
    LaserImpulse(
        double init_time,
        double init_phase,
        double intensity, 
        double freq, 
        double duration, 
        double step
        );
    ~LaserImpulse();
    SparseRMatrix getPMatr();
    SparseRMatrix getNMatr();
    SparseCMatrix& getImpulse(double t);
};

#endif