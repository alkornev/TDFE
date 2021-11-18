rm *.so *.o *.pyf
f2py --fcompiler=gfortran -m HermitSplines -c HermitSplines.f90
