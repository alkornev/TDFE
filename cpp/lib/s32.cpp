#include "s32.h"

#include <stdexcept>
#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>

/*
 * This module contains basis functions for S32 space for non-uniform grids
 */


void check_order(double x1i, double xi, double xi1){
    if (x1i > xi || xi > xi1 || x1i > xi1) {
        throw std::invalid_argument("x1i, xi, xi1 should be in increasing order.");
    }
}

CubicHermiteSplines::CubicHermiteSplines(const Eigen::Ref<const Eigen::VectorXd>& initGrid)
{
    grid = initGrid;
    nPoints = initGrid.size();
    splinesPerNode = 2;
    spaceDim = nPoints*splinesPerNode;
    initRefPoints();
}

CubicHermiteSplines::~CubicHermiteSplines(){}

void CubicHermiteSplines::initCollocGrid(int leftBC, int rightBC)
{
    double dt;

    double sh21, sh22;
    sh21 = (1.0 - 1.0/sqrt(3.0))*0.5;
    sh22 = (1.0 + 1.0/sqrt(3.0))*0.5;

    double sh31, sh32, sh33;
    sh31 = (1.0 - sqrt(0.6))*0.5;
    sh32 = 0.5;
    sh33 = (1.0 + sqrt(0.6))*0.5;

    int nCollocPoints = splinesPerNode * (nPoints-1);
    if (leftBC == -1) nCollocPoints += 1;
    if (rightBC == -1) nCollocPoints += 1;

    collocGrid = Eigen::VectorXd(nCollocPoints);
    dt = grid[1] - grid[0];
    int shift = 0;
    if (leftBC == -1){
        collocGrid[0] = grid[0] + dt * sh31;
        collocGrid[1] = grid[0] + dt * sh32;
        collocGrid[2] = grid[0] + dt * sh33;
        shift++;
    } else {
        collocGrid[0] = grid[0] + dt * sh21;
        collocGrid[1] = grid[0] + dt * sh22;
    }
    
    for(int i=1; i < nPoints-1; i++){
        dt = grid[i+1] - grid[i];
        collocGrid[2*i+shift] = grid[i] + dt * sh21;
        collocGrid[2*i+1+shift] = grid[i] + dt * sh22;
    }

    dt = grid[nPoints - 1] - grid[nPoints - 2];
    if (rightBC == -1){
        collocGrid[nCollocPoints-3] = grid[nPoints-2] + dt * sh31;
        collocGrid[nCollocPoints-2] = grid[nPoints-2] + dt * sh32;
        collocGrid[nCollocPoints-1] = grid[nPoints-2] + dt * sh33;
    } else {
        collocGrid[nCollocPoints-2] = grid[nPoints-2] + dt * sh21;
        collocGrid[nCollocPoints-1] = grid[nPoints-2] + dt * sh22;
    }
}

void CubicHermiteSplines::initRefPoints()
{
    int rangeNumber, splineNumber;
    double x1i, xi, xi1, res=0.0;

    leftPoints = Eigen::VectorXd(spaceDim);
    rightPoints = Eigen::VectorXd(spaceDim);
    midPoints = Eigen::VectorXd(spaceDim);

    for (int i = 0; i < spaceDim; i++)
    {
        rangeNumber = i / splinesPerNode;
        splineNumber = i % splinesPerNode;
        xi=grid[rangeNumber];
        if (rangeNumber != 0) x1i = grid[rangeNumber-1];
        else x1i = xi;

        if (rangeNumber < nPoints-1) xi1 = grid[rangeNumber+1];
        else xi1 = xi;
        
        leftPoints[i] = x1i;
        midPoints[i] = xi;
        rightPoints[i] = xi1;
    }
}

double CubicHermiteSplines::fBSpline(double t, int i)
{
    int rangeNumber, splineNumber;
    double x1i, xi, xi1, res = 0.0;
    rangeNumber = i / splinesPerNode;
    splineNumber = i % splinesPerNode;

    x1i = leftPoints[i];
    xi = midPoints[i];
    xi1 = rightPoints[i];


    switch (splineNumber)
    { 
        case 0: res=phi(t,x1i,xi,xi1); break;
        case 1: res=psi(t,x1i,xi,xi1); break;
    }

    return res;
}

double CubicHermiteSplines::d1BSpline(double t, int i)
{
    int rangeNumber, splineNumber;
    double x1i, xi, xi1, res = 0.0;
    rangeNumber = i / splinesPerNode;
    splineNumber = i % splinesPerNode;

    x1i = leftPoints[i];
    xi = midPoints[i];
    xi1 = rightPoints[i];

    switch (splineNumber)
    { 
        case 0: res=phi1(t,x1i,xi,xi1); break;
        case 1: res=psi1(t,x1i,xi,xi1); break;
    }
    return res;
}

double CubicHermiteSplines::d2BSpline(double t, int i)
{
    int rangeNumber, splineNumber;
    double x1i, xi, xi1, res = 0.0;
    rangeNumber = i / splinesPerNode;
    splineNumber = i % splinesPerNode;

    x1i = leftPoints[i];
    xi = midPoints[i];
    xi1 = rightPoints[i];

    switch (splineNumber)
    { 
        case 0: res=phi2(t,x1i,xi,xi1); break;
        case 1: res=psi2(t,x1i,xi,xi1); break;
    }
    return res;
}

double CubicHermiteSplines::phi(double x, double x1i, double xi, double xi1)
{
    check_order(x1i, xi, xi1);
    
    double hi;
    double lp;
    double rp;
    if (x1i <= x && x < xi) {
        lp = x - x1i;
        rp = x - xi;
        hi = xi - x1i;
        return -2.0*lp*lp*(rp-0.5*hi)/(hi*hi*hi);
    } else if (xi < x && x <= xi1) {
        lp = x - xi;
        rp = x - xi1;
        hi = xi1 - xi;
        return 2.0*rp*rp*(lp+0.5*hi)/(hi*hi*hi);
    } else if (xi == x) {
        return 1.0;
    } else {
        return 0.0;
    }
}

double CubicHermiteSplines::phi1(double x, double x1i, double xi, double xi1)
{
    check_order(x1i, xi, xi1);
    
    double hi;
    double lp;
    double rp;
    if (x1i <= x && x < xi) {
        lp = x - x1i;
        rp = x - xi;
        hi = xi - x1i;
        return -2.0*(2*lp*rp-lp*hi+lp*lp)/(hi*hi*hi);
    } else if (xi < x && x <= xi1) {
        lp = x - xi;
        rp = x - xi1;
        hi = xi1 - xi;
        // 2 * (DR * DR + (DL + 0.5 * HK) * DR + (DL + 0.5 * HK) * DR) / (HK * HK * HK);
        return 2.0*(2*rp*lp+rp*hi+rp*rp)/(hi*hi*hi);
    } else {
        return 0.0f;
    }
}

double CubicHermiteSplines::phi2(double x, double x1i, double xi, double xi1)
{
    check_order(x1i, xi, xi1);   
    

    double hi;
    double lp;
    double rp;
    if (x1i <= x && x < xi) {
        lp = x - x1i;
        rp = x - xi;
        hi = xi - x1i;
        return -4.0*(rp-0.5*hi+2.0*lp)/(hi*hi*hi);
    } else if (xi < x && x <= xi1) {
        lp = x - xi;
        rp = x - xi1;
        hi = xi1 - xi;
        return 4.0*(lp+0.5*hi+2.0*rp)/(hi*hi*hi);
    } else {
        return 0.0f;
    }
}

double CubicHermiteSplines::psi(double x, double x1i, double xi, double xi1)
{ 
    check_order(x1i, xi, xi1);
    
    double hi;
    double lp;
    double rp;
    if (x1i <= x && x < xi) {
        lp = x - x1i;
        rp = x - xi;
        hi = xi - x1i;
        return lp*lp*rp / (hi*hi);
    } else if (xi < x && x <= xi1) {
        lp = x - xi;
        rp = x - xi1;
        hi = xi1 - xi;
        return rp*rp*lp / (hi*hi);
    } else {
        return 0.0f;
    }
}

double CubicHermiteSplines::psi1(double x, double x1i, double xi, double xi1)
{
    check_order(x1i, xi, xi1);
    
    double hi;
    double lp;
    double rp;
    if (x1i <= x && x < xi) {
        lp = x - x1i;
        rp = x - xi;
        hi = xi - x1i;
        return (2*lp*rp + lp*lp) / (hi*hi);
    } else if (xi < x && x <= xi1) {
        lp = x - xi;
        rp = x - xi1;
        hi = xi1 - xi;
        return (2*rp*lp + rp*rp) / (hi*hi);
    } else if (xi == x){
        return 1.0f;
    } else {
        return 0.0f;
    }
}

double CubicHermiteSplines::psi2(double x, double x1i, double xi, double xi1)
{
    check_order(x1i, xi, xi1);

    double lp;
    double rp;
    double hi;
    if (x1i <= x && x < xi) {
        hi = xi - x1i;
        lp = x - x1i;
        rp = x - xi;
        return 2.0f*(2.0f*lp + rp) / (hi*hi);
    } else if (xi < x && x <= xi1) {
        hi = xi1 - xi;
        lp = x - xi;
        rp = x - xi1;
        return 2.0f*(2.0f*rp + lp) / (hi*hi);
    } else {
        return 0.0f;
    }
}


CHermiteBC::CHermiteBC(const Eigen::Ref<const Eigen::VectorXd>& initGrid, int lBC, int rBC) : space(initGrid)
{
    space = CubicHermiteSplines(initGrid);
    leftBC = lBC;
    rightBC = rBC;
    //std::cout << space.spaceDim << std::endl;
    splineBCdim = space.spaceDim - (int) (leftBC >= 0) - (int) (rightBC >= 0);
    
    //std::cout << splineBCdim << std::endl;
    space.initCollocGrid(leftBC, rightBC);
    initOverlapMatrix();
}

CHermiteBC::~CHermiteBC(){}

Eigen::SparseMatrix<double> CHermiteBC::generateSPMatr(){
    int dim = splineBCdim;
    int shift = (int) (leftBC >= 0);

    Eigen::SparseMatrix<double> pMatrix(dim, dim);

    for (int i = 0; i < splineBCdim; i++){
        double xi = space.collocGrid[i];
        for (int j = 0; j < splineBCdim; j++){  
            double leftJ = space.leftPoints[j+shift];
            double rightJ = space.rightPoints[j+shift];
            if (fBSplineBC(xi, j) != 0.0) {
                pMatrix.insert(i, j) = fBSplineBC(xi, j);
            }
        }
    }

    return pMatrix;
}

Eigen::SparseMatrix<double> CHermiteBC::generateSNMatr(){
    int dim = splineBCdim;
    int shift = (int) (leftBC >= 0);

    Eigen::SparseMatrix<double> nMatrix(dim, dim);

    for (int i = 0; i < splineBCdim; i++){
        // std::cout << "axi:" << axi << " bxi:" << bxi << std::endl;
        for (int j = 0; j < splineBCdim; j++){   
            if (overlapMatrix(i, j) != 0.0) {
                nMatrix.coeffRef(i, j) = overlapMatrix(i,j);
            }
        }
    }

    return nMatrix;
}


void CHermiteBC::initOverlapMatrix()
{
    overlapMatrix = Eigen::MatrixXd::Zero(splineBCdim, splineBCdim);
    
    double gaussPoints[3]={0.2386191860831969086305017,0.6612093864662645136613996,0.9324695142031520278123016};
    double gaussWeights[3]={0.4679139345726910473898703,0.3607615730481386075698335,0.1713244923791703450402961};
    
    double leftI, leftJ, rightI, rightJ, midI, midJ;
    double left, right, middle;
    double a, b, A, B, t;
    int shift = (int) (leftBC >= 0);
    
    for (int i=0; i < overlapMatrix.rows(); i++)
    {
        leftI=space.leftPoints[i+shift];
        midI=space.midPoints[i+shift];
        rightI=space.rightPoints[i+shift];
        for (int j=0; j < overlapMatrix.cols(); j++)
        {
            leftJ=space.leftPoints[j+shift];
            midJ=space.midPoints[j+shift];
            rightJ=space.rightPoints[j+shift];
            
            overlapMatrix(i, j)=0.0;
            if ((leftI>rightJ) || (leftJ>rightI))
                overlapMatrix(i, j)=0.0;
            else
            {
                left=std::max(leftI,leftJ);
                right=std::min(rightI,rightJ);
                if ((leftI==leftJ) && (rightI==rightJ))
                {
                    left=leftI;
                    middle=space.midPoints[i+shift];
                    right=rightI;
                    overlapMatrix(i, j)=0.0;
                    a=left;
                    b=middle;
                    A=0.5*(b-a);
                    B=0.5*(b+a);
                    for (int k=0; k<3; k++)
                    {
                        t=B+A*gaussPoints[k];
                        overlapMatrix(i,j)+=A*gaussWeights[k]*fBSplineBC(t, i)*fBSplineBC(t, j);
                        t=B-A*gaussPoints[k];
                        overlapMatrix(i,j)+=A*gaussWeights[k]*fBSplineBC(t, i)*fBSplineBC(t, j);
                    }
                    a=middle;
                    b=right;
                    A=0.5*(b-a);
                    B=0.5*(b+a);
                    for (int k=0; k<3; k++)
                    {
                        t=B+A*gaussPoints[k];
                        overlapMatrix(i,j)+=A*gaussWeights[k]*fBSplineBC(t, i)*fBSplineBC(t, j);
                        t=B-A*gaussPoints[k];
                        overlapMatrix(i,j)+=A*gaussWeights[k]*fBSplineBC(t, i)*fBSplineBC(t, j);
                    }
                }
                else
                {
                    a=left;
                    b=right;
                    A=0.5*(b-a);
                    B=0.5*(b+a);
                    for (int k=0; k<3; k++)
                    {
                        t=B+A*gaussPoints[k];
                        overlapMatrix(i,j)+=A*gaussWeights[k]*fBSplineBC(t, i)*fBSplineBC(t, j);
                        t=B-A*gaussPoints[k];
                        overlapMatrix(i,j)+=A*gaussWeights[k]*fBSplineBC(t, i)*fBSplineBC(t, j);
                    }
                }
            }
        }
    }
}

double CHermiteBC::fBSplineBC(double t, int i)
{
    double res = 0.0;
    int shift = leftBC >= 0;
    
    if (i == 0 && leftBC >= 0) {
        return (int)(leftBC == 0)*space.fBSpline(t, 0) + (int)(leftBC == 1)*space.fBSpline(t, 1);
    } 
    else if (i == splineBCdim-1 && rightBC >= 0) {
        return (int)(rightBC == 0)*space.fBSpline(t, space.spaceDim-2) + (int)(rightBC == 1)*space.fBSpline(t, space.spaceDim-1);
    } else 
        return space.fBSpline(t, i+shift);
}

double CHermiteBC::d1BSplineBC(double t, int i)
{
    double res = 0.0;
    int shift = leftBC >= 0;
    
    if (i == 0 && leftBC >= 0) {
        return (int)(leftBC == 0)*space.d1BSpline(t, 0) + (int)(leftBC == 1)*space.d1BSpline(t, 1);
    } 
    else if (i == splineBCdim-1 && rightBC >= 0) {
        return (int)(rightBC == 0)*space.d1BSpline(t, space.spaceDim-2) + (int)(rightBC == 1)*space.d1BSpline(t, space.spaceDim-1);
    } else 
        return space.d1BSpline(t, i+shift);
}

double CHermiteBC::d2BSplineBC(double t, int i)
{
    double res = 0.0;
    int shift = leftBC >= 0;
    
    if (i == 0 && leftBC >= 0) {
        return (int)(leftBC == 0)*space.d2BSpline(t, 0) + (int)(leftBC == 1)*space.d2BSpline(t, 1);
    } 
    else if (i == splineBCdim-1 && rightBC >= 0) {
        return (int)(rightBC == 0)*space.d2BSpline(t, space.spaceDim-2) + (int)(rightBC == 1)*space.d2BSpline(t, space.spaceDim-1);
    } else 
        return space.d2BSpline(t, i+shift);
}

int CHermiteBC::getSpaceDim(){ return splineBCdim; }

Eigen::VectorXd CHermiteBC::getLeftPoints() { return space.leftPoints; }
Eigen::VectorXd CHermiteBC::getRightPoints() { return space.rightPoints; }
Eigen::VectorXd CHermiteBC::getMidPoints() { return space.midPoints; }
Eigen::VectorXd CHermiteBC::getCollocPoints() { return space.collocGrid; }

double CHermiteBC::locate(const Eigen::Ref<const Eigen::VectorXd>& coefs, double x)
{
    Eigen::VectorXd& grid = space.grid;
    int right = grid.size() - 1;
    int left = 0;
    int mid;

    if (x > grid[right] || x < grid[left]){
         throw std::invalid_argument("argument outside boundaries");
    }

    while (right - left > 1){
        mid = left + (right - left) / 2;
        if (grid[mid] < x) {
            left = mid;
        } else if (grid[mid] > x) {
            right = mid;
        } else {
            break;
        }
    }

    double result = 0.0;
    for(int i = left; left < std::min(left + 4, splineBCdim); i++){
        result += coefs[i]*fBSplineBC(x, i);
    }
            
    return result;
}
