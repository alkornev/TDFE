#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif
#include "hamiltonian3d.h"
#include "s32.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <Eigen/QR>

#include <cmath>
#include <array>
#include <iostream>
#include <complex>
#include <typeinfo>


Hamiltonian3D::Hamiltonian3D(const Eigen::Ref<const Eigen::VectorXd>& aGrid, 
    const Eigen::Ref<const Eigen::VectorXd>& bGrid, 
    const std::array<int, 4>& BCs,
    const std::array<double, 3>& masses, double regularization, int n) : aSplines(aGrid, BCs[0], BCs[1]), bSplines(bGrid, BCs[2], BCs[3])
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
        
        //std::cout << "Initialization..." << std::endl;
        nMatr = generateNMatr();
        //std::cout << "Normalization Matrix is ready! Continuing..." << std::endl;

        pMatr = generatePMatr();

        //std::cout << "Collocation Matrix is ready! Continuing..." << std::endl;

        pMatrInv = pMatr.inverse();

        //std::cout << "Inversion of Collocation Matrix is ready! Continuing..." << std::endl;

        h = generateTheHamiltonian();
        //std::cout << "Hamiltonian is ready! Continuing..." << std::endl;
        //getTheSpectrum();  
    }
Hamiltonian3D::~Hamiltonian3D(){}

Eigen::MatrixXd Hamiltonian3D::generateNMatr(){
    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;
    Eigen::MatrixXd nMatrix(nMax, nMax);

    for (int i = 0; i < nMax; i++){
        for (int j = 0; j < nMax; j++){   
            nMatrix(i, j) = aSplines.overlapMatrix(i / bNMax, j / bNMax) * bSplines.overlapMatrix(i % bNMax, j % bNMax);
        }
    }
    return nMatrix; 
}

Eigen::MatrixXd Hamiltonian3D::generatePMatr(){
    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;
    Eigen::MatrixXd pMatrix(nMax, nMax);

    for (int i = 0; i < nMax; i++){
        double axi = aSplines.space.collocGrid[i / bNMax];
        double bxi = bSplines.space.collocGrid[i % bNMax];
        // std::cout << "axi:" << axi << " bxi:" << bxi << std::endl;
        for (int j = 0; j < nMax; j++){   
            pMatrix(i, j) = aSplines.fBSplineBC(axi, j / bNMax) * bSplines.fBSplineBC(bxi, j % bNMax);
        }
    }
    return pMatrix;   
}

Eigen::MatrixXd Hamiltonian3D::generateTheHamiltonian(){
    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;
    Eigen::MatrixXd ham(nMax, nMax);
    double potential = 0.0;
    double laplacian = 0.0;
    double mu = 2.0*m[0]*m[1] / (m[0] + m[1]);
    double mu123 = 2.0*((m[0]+m[1])*m[2])/(m[0] + m[1] + m[2]);
    for (int i = 0; i < nMax; i++){
        double axi = aSplines.space.collocGrid[i / bNMax];
        double bxi = bSplines.space.collocGrid[i % bNMax];
        for (int j = 0; j < nMax; j++){   
            laplacian = -1.0/mu*aSplines.d2BSplineBC(axi, j / bNMax)*bSplines.fBSplineBC(bxi, j % bNMax) 
                        -1.0/mu123*aSplines.fBSplineBC(axi, j / bNMax)*bSplines.d2BSplineBC(bxi, j % bNMax);
            potential = signs[0]/sqrt(axi*axi + r) +
                        signs[1]/sqrt((bxi - m[0]*axi/(m[0]+m[1]))*(bxi - m[0]*axi/(m[0]+m[1])) + r) +
                        signs[2]/sqrt((bxi + m[1]*axi/(m[0]+m[1]))*(bxi + m[1]*axi/(m[0]+m[1])) + r);
            ham(i, j) = laplacian + 0.5*potential*aSplines.fBSplineBC(axi, j / bNMax)*bSplines.fBSplineBC(bxi, j % bNMax);
        }
    }
    ham = pMatrInv*ham;
    return ham;   
}

double Hamiltonian3D::getEigenfunction(const Eigen::Ref<const Eigen::VectorXd>& coefs, double x, double y)
{
    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;

    Eigen::VectorXd t = Eigen::VectorXd(nMax);
    //std::cout <<"size of points: "<< t.size() << std::endl;
    for(int j=0; j < nMax; j++){
        t[j] = aSplines.fBSplineBC(x, j / bNMax)*bSplines.fBSplineBC(y, j % bNMax);
    }


    //std::cout <<"size of coefs" << cshaped.size() << std::endl;

    double res = t.dot(coefs);

    //std::cout << res << std::endl;

    return res;
}

groundState Hamiltonian3D::getExpGroundState(double err, double dt){
    // calc exponent using pade approximation
    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;

    Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(nMax, nMax);
    Eigen::VectorXd eigVector = Eigen::VectorXd::Random(nMax);

    Eigen::MatrixXd S = eye - 0.5 * dt * h;
    Eigen::MatrixXd F = eye + 0.5 * dt * h;
    Eigen::PartialPivLU<Eigen::MatrixXd> dec(F);
    
    double eigValue = eigVector.dot(h*eigVector)/(eigVector.dot(eigVector));
    double residue = (h*eigVector - eigValue*eigVector).norm(); 

    while (residue > err){
        Eigen::VectorXd y = dec.solve(eigVector);
        eigVector = S*y;
        eigValue = eigVector.dot(h*eigVector)/eigVector.dot(eigVector);
        eigVector = eigVector/(eigVector.norm());
        residue = (h*eigVector - eigValue*eigVector).norm(); 
        //std::cout << eigValue << "    " << residue << std::endl;
    }
    struct groundState gs = {eigValue, eigVector};

    return gs;
}

void Hamiltonian3D::getTheSpectrum(){
    Eigen::EigenSolver<Eigen::MatrixXd> es;
  
    //int const n = Eigen::nbThreads( );
    //std::cout << "#Threads: " << n << std::endl;

    es.compute(h);
    
    this->eigVectors = es.eigenvectors().real();
    this->eigVals = es.eigenvalues();

    Eigen::VectorXi indices(eigVals.size());
    for (int i = 0 ; i != indices.size() ; i++) {
        indices[i] = i;
    };

    std::sort(indices.begin(), indices.end(), [&](const int& a, const int& b)->bool {return (eigVals[a].real() < eigVals[b].real());});

    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> P;
    P.indices() = indices;
    eigVectors = eigVectors * P;
    eigVals = eigVals.transpose() * P;


    for(int i=0; i<eigVals.size(); i++){
        double norm = sqrt(eigVectors.col(i).dot(nMatr*eigVectors.col(i)));
        eigVectors.col(i) /=norm;
        double phasenorm = getEigenfunction(this->eigVectors.col(i), 1.0, 1.0);
        phasenorm = phasenorm/abs(phasenorm);
        eigVectors.col(i) /= phasenorm;
    }
}

Eigen::VectorXcd Hamiltonian3D::getEigenvalues(){ return eigVals; }
Eigen::MatrixXd Hamiltonian3D::getEigenvectors(){ return eigVectors; }

Eigen::MatrixXd Hamiltonian3D::getPMatr(){ return pMatr; }
Eigen::MatrixXd Hamiltonian3D::getNMatr(){ return nMatr; }

Eigen::MatrixXd Hamiltonian3D::getHamiltonian(){ return h; }

Eigen::MatrixXd Hamiltonian3D::getTDHamiltonian(double t, double V){
    double w0 = 0.058;
    double I0 = 350;
    double T = 248;
    double V0 = V;
    double envArg = -t/T;
    int aNMax = aSplines.splineBCdim;
    int bNMax = bSplines.splineBCdim;
    int nMax = aNMax * bNMax;
    Eigen::MatrixXd electricField = Eigen::MatrixXd(nMax, nMax);
    double bxi = 0.0;
    double axi = 0.0;

    double pulse = V0*std::exp(-envArg*envArg)*sin(w0*t);
    double basisElement = 0.0;
    for(int i=0; i < nMax; i++){
        axi = aSplines.space.collocGrid[i / bNMax];
        bxi = bSplines.space.collocGrid[i % bNMax];
        for (int j=0; j < nMax; j++){
            basisElement = aSplines.fBSplineBC(axi, i / bNMax)*bSplines.fBSplineBC(bxi, i % bNMax);
            electricField(i, j) = (1.0+m[2]/(m[0]+m[1]+m[2]))*bxi*pulse*basisElement;
        }
    }
    return pMatrInv*electricField;
}
Eigen::VectorXcd Hamiltonian3D::evolutionStep(Eigen::VectorXcd state, int iter, double dt, double V){
    /*
        We want to calculate this formula:
        (I - i/2*A)*inv(I + i/2*A) x = c
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


    std::complex<double> im(0.0, 1.0);
    Eigen::MatrixXcd eye = Eigen::MatrixXcd::Identity(nMax, nMax);
    
    Eigen::MatrixXcd S = eye - 0.25 * dt * im * h.cast<std::complex<double>>();
    Eigen::MatrixXcd F = eye + 0.25 * dt * im * h.cast<std::complex<double>>();

    Eigen::PartialPivLU<Eigen::MatrixXcd> dec(F);

    Eigen::VectorXcd y;
    Eigen::MatrixXd idtV;
    Eigen::MatrixXcd SV;
    Eigen::MatrixXcd FV;

    y = dec.solve(state);
    y = S*y;

    Eigen::MatrixXcd Vt = getTDHamiltonian(0.5*iter*dt, V);

    SV = eye - 0.5 * dt * im * Vt.cast<std::complex<double>>();
    FV = eye + 0.5 * dt * im * Vt.cast<std::complex<double>>();
    
    Eigen::PartialPivLU<Eigen::MatrixXcd> decV(FV);
    y = decV.solve(y);
    y = SV*y;

    y = dec.solve(y);
    y = S*y;
    
    
    return y;
}

