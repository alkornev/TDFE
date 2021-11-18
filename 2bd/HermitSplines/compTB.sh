rm tbSplines.so *.o *.pyf
f2py --fcompiler=gfortran -m tbSplines -c tbSplines.f90
