#include "hamiltonian2d.h"
#include "s32.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>


Hamiltonian2D::Hamiltonian2D(const Eigen::Ref<const Eigen::VectorXd>& initGrid, 
                            int leftBC, int rightBC, 
                            const std::array<double, 2>& masses, double reg) : spl(initGrid, leftBC, rightBC) {
    m = masses;
    r = reg;
    nMatr = spl.overlapMatrix;
    pMatr = generatePMatr();
    potential = generatePotential();
    d2Matr = generateD2Matr();
    pMatrInv = pMatr.inverse();
    h2 = generateTheHamiltonian();
    //getTheSpectrum();
}

Hamiltonian2D::~Hamiltonian2D(){}

Eigen::MatrixXd Hamiltonian2D::generatePMatr() {
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

Eigen::MatrixXd Hamiltonian2D::getPMatr(){
    return pMatr;
}

Eigen::MatrixXd Hamiltonian2D::getNMatr(){
    return nMatr;
}

Eigen::MatrixXd Hamiltonian2D::getD2Matr(){
    return d2Matr;
}

Eigen::MatrixXd Hamiltonian2D::getPotential(){
    return potential;
}

Eigen::MatrixXd Hamiltonian2D::generateD2Matr(){
    int nMax = spl.splineBCdim;
    Eigen::MatrixXd d2(nMax, nMax);
    for(int i = 0; i < nMax; i++){
        double xi = spl.space.collocGrid[i];
        for(int j=0; j < nMax; j++){
            d2(i,j)=spl.d2BSplineBC(xi, j);
        }
    }

    return d2;
}

Eigen::MatrixXd Hamiltonian2D::generatePotential(){
    int nMax = spl.splineBCdim;
    Eigen::MatrixXd pot(nMax, nMax);
    for(int i = 0; i < nMax; i++){
        double xi = spl.space.collocGrid[i];
        for(int j=0; j < nMax; j++){
            pot(i,j)= -0.5/sqrt(xi*xi+r)*spl.fBSplineBC(xi, j);
        }
    }

    return pot;
}

Eigen::MatrixXd Hamiltonian2D::generateTheHamiltonian() {
    int nMax = spl.splineBCdim;
    Eigen::MatrixXd ham(nMax, nMax);
    double mu = m[0]*m[1]/(m[0] + m[1]);
    for(int i = 0; i < nMax; i++){
        double xi = spl.space.collocGrid[i];
        for(int j=0; j < nMax; j++){
            ham(i,j)=-1.0/(2.0*mu)*spl.d2BSplineBC(xi, j) - 0.5/sqrt(xi*xi+r)*spl.fBSplineBC(xi, j);
        }
    }
    ham = pMatrInv*ham;
    return ham;
}

double Hamiltonian2D::getEigenfunction(const Eigen::Ref<const Eigen::VectorXd>& coefs, double x) {
    Eigen::VectorXd cshaped = coefs.real();
    Eigen::VectorXd t = Eigen::VectorXd(spl.splineBCdim);
    for (int i=0; i < spl.splineBCdim; i++){
        t[i] = spl.fBSplineBC(x, i);
    }
    double res = t.dot(cshaped);
    return res;
}

void Hamiltonian2D::getTheSpectrum() {
    Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(h2);
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
Eigen::VectorXd Hamiltonian2D::getEigenvalues() { return eigVals; }
Eigen::MatrixXd Hamiltonian2D::getEigenvectors() { return eigVectors; }
