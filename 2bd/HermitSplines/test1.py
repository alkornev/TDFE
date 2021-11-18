#!/usr/bin/python
__author__="roudnev"
__date__ ="$Apr 9, 2009 2:49:25 PM$"
if __name__ == "__main__":
    from pylab import *
    from scipy import arange
    import numpy as np
    from HermitSplines import hermitsplines
#    print HermitSplines.__doc__
#subroutine InitDiscretization(NPoints,leftBound,rightBound,specialPoint, &
#      &                         densityLeft,densityRight,iBCLeft,iBCRight)

    y1i=-1;
    yi1=1;
    yi=0;
    grd=arange(y1i,yi1,0.0001)
    grid
    y=[hermitsplines.phi(x,y1i,yi,yi1) for x in grd]
    z=[hermitsplines.psi(x,y1i,yi,yi1) for x in grd]
    #w=[hermitsplines.phi1(x,y1i,yi,yi1) for x in grd]
    ##f=[hermitsplines.psi1(x,y1i,yi,yi1) for x in grd]
    #s=[hermitsplines.psi2(x,y1i,yi,yi1) for x in grd]
    #k=[hermitsplines.phi2(x,y1i,yi,yi1) for x in grd]
    #z=[hermitsplines.phi1(x,y1i,yi,yi1) for x in grd]
    plot(grd,y,grd,z)#,grd,w,grd,f,grd,s,grd,k)
    show()
