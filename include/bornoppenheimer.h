#ifndef BORNOPPENHEIMER2D_H
#define BORNOPPENHEIMER2D_H
#include "s32.h"
#include <Eigen/Dense>


class BornOppenheimer2D
{
public:
    CHermiteBC spl;
    std::array<double, 3> m;
    std::array<double, 3> signs;
    double r;
    double R;
    Eigen::MatrixXd h;
    Eigen::MatrixXd pMatr;
    Eigen::MatrixXd nMatr;
    Eigen::MatrixXd pMatrInv;
    Eigen::MatrixXd eigVectors;
    Eigen::VectorXd eigVals;
    BornOppenheimer2D(const Eigen::Ref<const Eigen::VectorXd>& initGrid, 
                                                            int leftBC, 
                                                            int rightBC, 
                                                            const std::array<double, 3>& masses,
                                                            double reg,
                                                            double nuclRad);
    ~BornOppenheimer2D();
    Eigen::MatrixXd generatePMatr();
    Eigen::MatrixXd generateTheHamiltonian();
    double getEigenfunction(const Eigen::Ref<const Eigen::VectorXd>& coefs, double x);
    void getTheSpectrum();
    Eigen::VectorXd getEigenvalues();
    Eigen::MatrixXd getEigenvectors();
    Eigen::MatrixXd getPMatr();
    Eigen::MatrixXd getNMatr();
};
#endif