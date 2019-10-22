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

sp_list = []
for i in range(nx):
    for j in range(ny):
        cc = tricubicMat[i,j]
        if cc == 0: continue
        if abs(cc) == 1: continue
        if float(abs(cc)) not in sp_list:
            sp_list.append(float(abs(cc)))
print('Number of unique constants: %d'%len(sp_list))
sp_list.sort()
c_tri_list = str(sp_list).replace('[','{').replace(']','}')
ff.write('    const double tri_consts[%d] = '%len(sp_list))
ff.write(c_tri_list)
ff.write(';\n')

ss = ' '
for i in range(nx):
#for i in range(1):
    flag = True
    for j in range(ny):
        cc = tricubicMat[i,j]
        
        if cc == 0:
            continue

        if flag:
            ss = ' '
            ff.write('    ')
        else:
            ss = '+'
        
        if abs(cc) == 1:
#            print('coefs[%d] %s= b[%d]; '%(i,ss,j),end='')
            if ss == ' ':
                if cc == 1:
                    ff.write('coefs[%d] = b[%d]; '%(i,j))
                else:
                    ff.write('coefs[%d] = -b[%d]; '%(i,j))
            else:
                if cc == 1:
                    ss = '+'
                else:
                    ss = '-'
                ff.write('coefs[%d] %s= b[%d]; '%(i,ss,j))
        else:
            ind = sp_list.index(float(abs(cc)))
            if ss == ' ':
                if cc > 0:
                    ff.write('coefs[%d] = tri_consts[%d] * b[%d]; '%(i,ind,j))
                else:
                    ff.write('coefs[%d] = -(tri_consts[%d] * b[%d]); '%(i,ind,j))
            else:
                if cc > 0:
                    ss = '+'
                else:
                    ss = '-'
                #ff.write('coefs[%d] %s= b[%d]; '%(i,ss,j))
                ff.write('coefs[%d] %s= tri_consts[%d] * b[%d]; '%(i,ss,ind,j))
        flag = False
#    print()
    ff.write('\n')

ff.write(
'''
    return coefs;
}
''')

