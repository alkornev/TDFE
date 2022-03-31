#include "bornoppenheimer.h"
#include "s32.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>


BornOppenheimer2D::BornOppenheimer2D(const Eigen::Ref<const Eigen::VectorXd>& initGrid, 
                            int leftBC, int rightBC, 
                            const std::array<double, 3>& masses, double reg, double nuclRad) : spl(initGrid, leftBC, rightBC) {
    r = reg;
    R = nuclRad;

    m = masses;

    nMatr = spl.overlapMatrix;
    pMatr = generatePMatr();
    pMatrInv = pMatr.inverse();
    h = generateTheHamiltonian();
    //getTheSpectrum();
}

BornOppenheimer2D::~BornOppenheimer2D(){}

Eigen::MatrixXd BornOppenheimer2D::generatePMatr() {
    int nMax = spl.splineBCdim;
    Eigen::MatrixXd pMatr(nMax, nMax);

    for (int i = 0; i < nMax; i++){
        double xi = spl.space.collocGrid[i];
        for (int j = 0; j < nMax; j++){   
            pMatr(i, j) = spl.fBSplineBC(xi, j);
        }
    }
    return pMatr;   
}

Eigen::MatrixXd BornOppenheimer2D::getPMatr(){
    return pMatr;
}

Eigen::MatrixXd BornOppenheimer2D::getNMatr(){
    return nMatr;
}

Eigen::MatrixXd BornOppenheimer2D::generateTheHamiltonian() {
    int nMax = spl.splineBCdim;
    Eigen::MatrixXd ham(nMax, nMax);
    double mu = 2*m[0]*m[1]/(m[0] + m[1]);
    double potential = 0.0;
    for(int i = 0; i < nMax; i++){
        double xi = spl.space.collocGrid[i];
        for(int j=0; j < nMax; j++){
            potential = -1.0/sqrt((xi - m[0]*R/(m[0]+m[1]))*(xi - m[0]*R/(m[0]+m[1])) + r)
                        -1.0/sqrt((xi + m[1]*R/(m[0]+m[1]))*(xi + m[1]*R/(m[0]+m[1])) + r);
            ham(i,j)=-1.0/2.0*spl.d2BSplineBC(xi, j) + 0.5*potential*spl.fBSplineBC(xi, j);
        }
    }
    ham = pMatrInv*ham;
    return ham;
}

double BornOppenheimer2D::getEigenfunction(const Eigen::Ref<const Eigen::VectorXd>& coefs, double x) {
    Eigen::VectorXd cshaped = coefs.real();
    Eigen::VectorXd t = Eigen::VectorXd(spl.splineBCdim);
    for (int i=0; i < spl.splineBCdim; i++){
        t[i] = spl.fBSplineBC(x, i);
    }
    double res = t.dot(cshaped);
    return res;
}

void BornOppenheimer2D::getTheSpectrum() {
    Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(h);
    this->eigVectors = es.eigenvectors().real();
    this->eigVals = es.eigenvalues().real();

    Eigen::VectorXi indices(eigVals.size());
    for (int i = 0 ; i != indices.size() ; i++) {
        indices[i] = i;
    };

    std::sort(indices.begin(), indices.end(), [&](const int& a, const int& b)->bool {return (eigVals[a] < eigVals[b]);});

    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> P;
    P.indices() = indices;
    eigVectors = eigVectors * P;
    eigVals = eigVals.transpose() * P;


    for(int i=0; i<eigVals.size(); i++){
        double norm = sqrt(eigVectors.col(i).dot(nMatr*eigVectors.col(i)));
        eigVectors.col(i) /=norm;
        double phasenorm = getEigenfunction(this->eigVectors.col(i), 1.0);
        phasenorm = phasenorm/abs(phasenorm);
        eigVectors.col(i) /= phasenorm;
    }
}
Eigen::VectorXd BornOppenheimer2D::getEigenvalues() { return eigVals; }
Eigen::MatrixXd BornOppenheimer2D::getEigenvectors() { return eigVectors; }
