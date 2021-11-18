#!/usr/bin/python
__author__="roudnev"
__date__ ="$Apr 9, 2009 2:49:25 PM$"
if __name__ == "__main__":
    from pylab import *
    from scipy import arange
    import numpy as np
    from HermitSplines import hermitsplines
#    print hermitsplines.__doc__
#subroutine InitDiscretization(NPoints,leftBound,rightBound,specialPoint, &
#      &                         densityLeft,densityRight,iBCLeft,iBCRight)
    grd=arange(-1.0,1.1,0.5)
    grd=np.sin(np.pi/2*grd)
    print "The generated Grid"
    print grd
    hermitsplines.initdiscretization( grd, -1,-1, 0)
    x=arange(-1.0,1.001,0.001)
    print hermitsplines.ncx
    print "The grid:"
    for t in hermitsplines.xg:
        print t,
    print "\nThe collocation points:"
    for t in hermitsplines.xgc:
        print t,
    print
    for i in xrange(hermitsplines.ncx):
        print hermitsplines.xgc[i],
    print
    for i in xrange(hermitsplines.ncx):
        print i,
        y=[hermitsplines.hsplinebc(i,t) for t in x]
        plot(x,y)
    z=np.zeros(hermitsplines.ncx)
    plot(hermitsplines.xgc,z,'.')
    show()
    
