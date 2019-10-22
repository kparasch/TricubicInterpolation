from __future__ import print_function
import numpy as np
import sys
sys.path.append('../pyTricubic')

import tricubic_matrix

tricubicMat = tricubic_matrix.tricubicMat
ff = open('coefficients_function.c','w')

ff.write(
'''#include "coefs.h"

double* tricubic_get_coefs(double* b)
{
    double* coefs = (double*)malloc(64*sizeof(double));

''')

nx = tricubicMat.shape[0]
ny = tricubicMat.shape[1]

ss = ' '
for i in range(nx):
    flag = True
    for j in range(ny):
        cc = tricubicMat[i,j]
        
        if cc == 0:
            continue

        if flag:
            ss = ' '
            ff.write('    ')
        else:
            if cc > 0:
                ss = '+'
            else:
                ss = '-'
        
        if abs(cc) == 1:
#            print('coefs[%d] %s= b[%d]; '%(i,ss,j),end='')
            ff.write('coefs[%d] %s= b[%d]; '%(i,ss,j))
        else:
#            print('coefs[%d] %s= %d * b[%d]; '%(i,ss,abs(cc),j),end='')
            ff.write('coefs[%d] %s= %d. * b[%d]; '%(i,ss,abs(cc),j))
        flag = False
#    print()
    ff.write('\n')

ff.write(
'''
    return coefs;
}
''')

