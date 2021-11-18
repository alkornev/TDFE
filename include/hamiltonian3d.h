#ifndef HAMILTONIAN3D_H
#define HAMILTONIAN3D_H

#include "s32.h"
#include "hamiltonian2d.h"
#include <Eigen/Dense>
#include <array>


struct groundState{
    double eigValue;
    Eigen::VectorXd eigVector;
};


class Hamiltonian3D
{
public:
    double r;
    CHermiteBC aSplines;
    CHermiteBC bSplines;
    Eigen::MatrixXd h;
    std::array<double, 3> m;  
    std::array<double, 3> signs;  
    Eigen::MatrixXd pMatr;
    Eigen::MatrixXd nMatr;
    Eigen::MatrixXd pMatrInv;
    Eigen::MatrixXd eigVectors;
    Eigen::VectorXcd eigVals;
    Hamiltonian3D(const Eigen::Ref<const Eigen::VectorXd>& aGrid, 
                    const Eigen::Ref<const Eigen::VectorXd>& bGrid, 
                    const std::array<int, 4>& BCs,
                    const std::array<double, 3>& masses,
                    double regularization, int n);
    ~Hamiltonian3D();
    Eigen::MatrixXd generateNMatr();
    Eigen::MatrixXd generatePMatr();
    Eigen::MatrixXd generateTheHamiltonian();
    double getEigenfunction(const Eigen::Ref<const Eigen::VectorXd>& coefs, double x, double y);
    void getTheSpectrum();
    Eigen::VectorXcd getEigenvalues();
    Eigen::MatrixXd getEigenvectors();
    Eigen::MatrixXd getPMatr();
    Eigen::MatrixXd getNMatr();
    Eigen::MatrixXd getHamiltonian();
    Eigen::MatrixXd getTDHamiltonian(double dt, double V);
    groundState getExpGroundState(double err, double dt);
    Eigen::VectorXcd evolutionStep(Eigen::VectorXcd state, int iters, double dt, double V);
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