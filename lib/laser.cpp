#include "laser.h"


SparseCMatrix& LaserImpulse::getImpulse(double t)
{
    double V0 = std::sqrt(I0/3.5095);
    double envArg = t/tau;

    int bNMax = bSplines.splineBCdim;
    SparseCMatrix ebField = SparseCMatrix(bNMax, bNMax);
    
    ebField.reserve(Eigen::VectorXi::Constant(bNMax,10));
    double bxi = 0.0;

    double pulse = V0*std::exp(-envArg*envArg)*cos(w0*(t + t0) + phi);
    
        //std::cout << V0 << " " << pulse << " " << w0 << " " << phi << std::endl;

    for(int i=0; i < bNMax; i++){
        bxi = bSplines.space.collocGrid[i];
        for (int j=0; j < bNMax; j++){
            if (bpMatr.coeff(i, j) != 0.0) {
                ebField.coeffRef(i, j) = -(1.0+m[2]/(m[0]+m[1]+m[2]))*bxi*pulse*bpMatr.coeff(i, j);
            }
        }
    }
    ebField.makeCompressed();

    SparseCMatrix eField = Eigen::kroneckerProduct(apMatr, ebField).eval();

    return eField;
}
