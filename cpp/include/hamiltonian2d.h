#ifndef HAMILTONIAN2D_H
#define HAMILTONIAN2D_H
#include "s32.h"
#include <Eigen/Dense>


class Hamiltonian2D
{
public:
    CHermiteBC spl;
    std::array<double, 2> m;
    double r;
    Eigen::MatrixXd h2;
    Eigen::MatrixXd pMatr;
    Eigen::MatrixXd nMatr;
    Eigen::MatrixXd d2Matr;
    Eigen::MatrixXd potential;
    Eigen::MatrixXd pMatrInv;
    Eigen::MatrixXd eigVectors;
    Eigen::VectorXd eigVals;
    Hamiltonian2D(const Eigen::Ref<const Eigen::VectorXd>& initGrid, int leftBC, int rightBC, const std::array<double, 2>& masses, double reg);
    ~Hamiltonian2D();
    Eigen::MatrixXd generatePMatr();
    Eigen::MatrixXd generateD2Matr();
    Eigen::MatrixXd generatePotential();
    Eigen::MatrixXd generateTheHamiltonian();
    double getEigenfunction(const Eigen::Ref<const Eigen::VectorXd>& coefs, double x);
    void getTheSpectrum();
    Eigen::VectorXd getEigenvalues();
    Eigen::MatrixXd getEigenvectors();
    Eigen::MatrixXd getPMatr();
    Eigen::MatrixXd getNMatr();
    Eigen::MatrixXd getD2Matr();
    Eigen::MatrixXd getPotential();
};
#endif