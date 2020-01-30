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
    double b_val;
    int ii;

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

scale = ['1.', 'dx', 'dy', 'dz', '(dx * dy)', '(dx * dz)', '(dy * dz)', '(dx * dy) * dz']
iistrings = ['    ii =  8 * ( (ix  ) + nx * ( (iy  ) + ny * (iz  ) ) );\n',
             '    ii =  8 * ( (ix+1) + nx * ( (iy  ) + ny * (iz  ) ) );\n',
             '    ii =  8 * ( (ix  ) + nx * ( (iy+1) + ny * (iz  ) ) );\n',
             '    ii =  8 * ( (ix+1) + nx * ( (iy+1) + ny * (iz  ) ) );\n',
             '    ii =  8 * ( (ix  ) + nx * ( (iy  ) + ny * (iz+1) ) );\n',
             '    ii =  8 * ( (ix+1) + nx * ( (iy  ) + ny * (iz+1) ) );\n',
             '    ii =  8 * ( (ix  ) + nx * ( (iy+1) + ny * (iz+1) ) );\n',
             '    ii =  8 * ( (ix+1) + nx * ( (iy+1) + ny * (iz+1) ) );\n']

cgroups = [ 
#           [42,43,46,47,58,59,62,63], #64
#           [0,1,4,5,16,17,20,21]  #1line
#           [2,3,6,7,8,9,10,11,12,13,14]   #g1 
#           [18,19,22,23,24,25,26,27,28,29,30,31] #g2
            [32,33,34,35,48,49,50,51] #g3
          ]
bgroups = []
for ii in cgroups[0]:
    for jj in range(64):
        if tricubicMat[ii,jj]:
            if jj not in bgroups:
                bgroups.append(jj)
print(bgroups)
flags = np.ones([nx])
flags = np.zeros([nx])
ss = ' '
for gg in range(7):
    if gg == 1: break
    for kk in range(8):
        ff.write(iistrings[kk])
        for ll in range(8):
            if ll > 0:
                ff.write('    ii++;\n')
            j = 8*ll + kk
            if j not in bgroups: 
                continue
            ff.write('    b_val = LUT[ii];\n')
            for i in range(nx):
                if i not in cgroups[gg]:
                    continue

                cc = tricubicMat[i,j]
                
                if cc == 0:
                    continue
    
                if flags[i]:
                    ss = ' '
    #                ff.write('    ')
                else:
                    ss = '+'
                
                if abs(cc) == 1:
    #                print('coefs[%d] %s= b[%d]; '%(i,ss,j),end='')
                    if ss == ' ':
                        if cc == 1:
                            ff.write('    coefs[%d] = b_val;\n'%(i))
                        else:
                            ff.write('    coefs[%d] = -b_val;\n'%(i))
                    else:
                        if cc == 1:
                            ss = '+'
                        else:
                            ss = '-'
                        ff.write('    coefs[%d] %s= b_val;\n'%(i,ss))
                else:
                    ind = sp_list.index(float(abs(cc)))
                    if ss == ' ':
                        if cc > 0:
                            ff.write('    coefs[%d] = tri_consts[%d] * b_val;\n'%(i,ind))
                        else:
                            ff.write('    coefs[%d] = -(tri_consts[%d] * b_val);\n'%(i,ind))
                    else:
                        if cc > 0:
                            ss = '+'
                        else:
                            ss = '-'
                        #ff.write('coefs[%d] %s= b[%d]; '%(i,ss,j))
                        ff.write('    coefs[%d] %s= tri_consts[%d] * b_val;\n'%(i,ss,ind))
                flags[i] = 0
    #        print()
    #        ff.write('\n')

ff.write(
'''
    return coefs;
}
''')

