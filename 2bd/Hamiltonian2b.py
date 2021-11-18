#! /usr/bin/python
__author__="roudnev"
__date__ ="$Apr 6, 2009 8:01:10 PM$"
"""
This module configures the parameters of discretization for the few-body system.
"""
from HermitSplines.HermitSplines import hermitsplines
import math
from scipy import arange, zeros,linspace
from scipy import interpolate
import scipy.special
import scipy.integrate
import scipy.linalg
import matplotlib.pyplot as plt

# Physical parameters
# Number of bound states to reproduce
NBound=1
class Hamiltonian2B:
    def __init__(self,grid,iBCLeft=0,iBCRight=1,cc=1.0):
        # Discretization
        self.Npoints=len(grid)-1
        hermitsplines.initdiscretization(grid,iBCLeft,iBCRight)
        self.Nmax=len(hermitsplines.xgc)
        self.colpoints=hermitsplines.xgc
        # Masses of the particles, a.u.
        m1=1.0
#        m1=4.002603
#        m2=4.002603
#        m1=19.9924
        m2=m1*1822.888485
#        self.mu=2*1822.888485*m1*m2/(m1+m2)
        self.mu=m1*m2/(m1+m2)
        self.nMatr=self.generateNMatr()
        self.pMatr=self.generatePMatr()
        self.pMatrInv=scipy.linalg.inv(self.pMatr)
        self.h2=self.generateTheHamiltonian()


    def generatePMatr(self):
        res=scipy.array( [[hermitsplines.hsplinebc(i,self.colpoints[j]) for i in xrange(self.Nmax)] for j in xrange(self.Nmax)],'d')#.transpose()
        return res

    def generateNMatr(self):
        res=scipy.zeros((self.Nmax,self.Nmax),'d')
        Pn=hermitsplines.hsplinebc
        for j in xrange(self.Nmax):
            leftJ=hermitsplines.leftbounds[j]
            rightJ=hermitsplines.rightbounds[j]
            for i in xrange(self.Nmax):
              leftI=hermitsplines.leftbounds[i]
              rightI=hermitsplines.rightbounds[i]
              if (leftI>rightJ) or (leftJ>rightI):
                  res[i,j]=0.0
              else:
                  left=max(leftI,leftJ)
                  right=min(rightI,rightJ)
                  res[i,j]=scipy.integrate.quad(lambda x:Pn(i,x)*Pn(j,x),left,right)[0]
        return res
    #
    #  Now we should construct the 2-body Hamiltonian
    #
    def generateTheHamiltonian(self,cc=1.0):
          v2=lambda t:-1.0/t #(math.sqrt(t**2+0.01))
          spl=hermitsplines.hsplinebc
          spl1=hermitsplines.hsplinebc1
          spl2=hermitsplines.hsplinebc2
          h=scipy.zeros((self.Nmax,self.Nmax),'d')
          for i in xrange(self.Nmax):
              xi=self.colpoints[i]
              for j in xrange(self.Nmax):
                  h[i,j]=-spl2(j,xi)/(2*self.mu)+cc*spl(j,xi)*v2(xi) #v2(xi)#*12.11928/12.11930#(12.11928/12.12)#*3.1668154e-6/3.1669e-6  #/0.5291772108)/12.12/3.1668154e-6
          h=scipy.dot(self.pMatrInv,h)
          return h

    def getDerivMatr(self):
          spl2=hermitsplines.hsplinebc2
          d=scipy.zeros((self.Nmax,self.Nmax),'d')
          for i in xrange(self.Nmax):
              xi=self.colpoints[i]
              for j in xrange(self.Nmax):
                  d[i,j]=spl2(j,xi)
          #d=scipy.dot(self.pMatrInv,h)
          return d
    def splinef(self,coef,x):
        cshaped=coef
        Pn=hermitsplines.hsplinebc
        t=scipy.array([Pn(n,x) for n in xrange(self.Nmax)])
        res=scipy.dot(cshaped,t)
        return res

    def getTheSpectrum(self):
    #
        eval,evec=scipy.linalg.eig(self.h2)
        evalabs=eval
        ind=scipy.argsort(evalabs)
        eval=eval[ind]
        evec=evec[:,ind]
        for i in xrange(len(evec)):
            norm=scipy.dot(evec[:,i],scipy.dot(self.nMatr,evec[:,i]))
            evec[:,i]=evec[:,i]/math.sqrt(norm)
            phasenorm=self.splinef(evec[:,i], 1.0)
            phasenorm=phasenorm/abs(phasenorm)
            evec[:,i]=evec[:,i]/phasenorm
    #        print i,norm
        self.eigenvalues,self.eigenvectors=eval,evec
        return eval,evec


###########


##############################################################################3
if __name__ == "__main__":
    leftEnd=0.0
    rightEnd=40.0
    for n in [50,75,100,125]:
      grid=linspace(0,1.0,n)
      grid=rightEnd*grid**2+leftEnd
      h=Hamiltonian2B(grid,0,1)
      eval,evec= h.getTheSpectrum()
      print eval[:3]

    t = linspace(0, 60, 200)
    x = scipy.array([h.splinef(evec[:,2], i) for i in t])
    plt.plot(t,x)
    plt.show()
#    n=150
#    grid=linspace(0,1.0,n)
#    grid=rightEnd*grid**2.5
#    h=Hamiltonian2B(grid,0,1)
#    eval,evec= h.getTheSpectrum()
#    print eval[:3]
