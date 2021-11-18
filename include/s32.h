#ifndef S32_H
#define S32_H
#include <Eigen/Dense>


void check_order(double x1i, double xi, double xi1);


class CubicHermiteSplines
{
public:
    int nPoints;
    int splinesPerNode;
    int spaceDim;
    
    Eigen::VectorXd grid;
    Eigen::VectorXd leftPoints;
    Eigen::VectorXd rightPoints;
    Eigen::VectorXd midPoints;
    Eigen::VectorXd collocGrid;
    
    CubicHermiteSplines(const Eigen::Ref<const Eigen::VectorXd>& initGrid);

    ~CubicHermiteSplines();

    void initCollocGrid(int leftBC, int rightBC);

    void initRefPoints();

    void initOverlapMatrix();

    double fBSpline(double t, int i);

    double d1BSpline(double t, int i);

    double d2BSpline(double t, int i);

    double phi(double x, double x1i, double xi, double xi1);

    double phi1(double x, double x1i, double xi, double xi1);

    double phi2(double x, double x1i, double xi, double xi1);

    double psi(double x, double x1i, double xi, double xi1);

    double psi1(double x, double x1i, double xi, double xi1);

    double psi2(double x, double x1i, double xi, double xi1);
};

class CHermiteBC
{
public:
    int leftBC;
    int rightBC;
    int splineBCdim;
    CubicHermiteSplines space;
    
    Eigen::MatrixXd overlapMatrix;

    CHermiteBC(const Eigen::Ref<const Eigen::VectorXd>& initGrid, int leftBC, int rightBC);
    ~CHermiteBC();

    void initOverlapMatrix();
    double fBSplineBC(double t, int i);

    double d1BSplineBC(double t, int i);

    double d2BSplineBC(double t, int i);
    int getSpaceDim();

    Eigen::VectorXd getLeftPoints();
    Eigen::VectorXd getRightPoints();
    Eigen::VectorXd getMidPoints();
};

#endif
