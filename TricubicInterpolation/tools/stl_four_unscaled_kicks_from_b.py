from __future__ import print_function
import numpy as np
import sys
sys.path.append('../pyTricubic')

import tricubic_matrix

tricubicMat = tricubic_matrix.tricubicMat
ff = open('coefficients.h','w')

ff.write(
'''#ifndef SIXTRACKLIB_COMMON_BE_TRICUB_COEFFICIENTS_C99_H__
#define SIXTRACKLIB_COMMON_BE_TRICUB_COEFFICIENTS_C99_H__

/*
#if !defined(SIXTRL_NO_SYSTEM_INCLUDES)
    #include <stddef.h>
    #include <math.h>
#endif */
/* !defined(SIXTRL_NO_SYSTEM_INCLUDES) */

#if !defined(SIXTRL_NO_INCLUDES)
    #include "sixtracklib/common/definitions.h"
    #include "sixtracklib/common/be_tricub/be_tricub.h"
#endif /* !defined(SIXTRL_NO_INCLUDES) */

#if !defined( _GPUCODE ) && defined( __cplusplus )
extern "C" {
#endif /* !defined( _GPUCODE ) && defined( __cplusplus ) */

SIXTRL_STATIC SIXTRL_FN void NS(tricub_unscaled_kicks_from_b1)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    NS(be_tricub_real_t) x1,
    NS(be_tricub_real_t) y1,
    NS(be_tricub_real_t) z1,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* kicks);

SIXTRL_STATIC SIXTRL_FN void NS(tricub_unscaled_kicks_from_b2)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    NS(be_tricub_real_t) x1,
    NS(be_tricub_real_t) y1,
    NS(be_tricub_real_t) z1,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* kicks);

SIXTRL_STATIC SIXTRL_FN void NS(tricub_unscaled_kicks_from_b3)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    NS(be_tricub_real_t) x1,
    NS(be_tricub_real_t) y1,
    NS(be_tricub_real_t) z1,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* kicks);

SIXTRL_STATIC SIXTRL_FN void NS(tricub_unscaled_kicks_from_b4)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    NS(be_tricub_real_t) x1,
    NS(be_tricub_real_t) y1,
    NS(be_tricub_real_t) z1,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* kicks);

#if !defined( _GPUCODE ) && defined( __cplusplus )
}
#endif /* !defined( _GPUCODE ) && defined( __cplusplus ) */

/* ************************************************************************* */
/* ************************************************************************* */

#if !defined( _GPUCODE ) && defined( __cplusplus )
extern "C" {
#endif /* !defined( _GPUCODE ) && defined( __cplusplus ) */

SIXTRL_INLINE void NS(tricub_unscaled_kicks_from_b1)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    NS(be_tricub_real_t) x1,
    NS(be_tricub_real_t) y1,
    NS(be_tricub_real_t) z1,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* kicks
    )
{
//    NS(be_tricub_real_t) const tri_consts[9] = {2.0, 3.0, 4.0, 6.0, 8.0, 9.0, 12.0, 18.0, 27.0};
    NS(be_tricub_real_t) coef;
    NS(be_tricub_real_t) akick;

    NS(be_tricub_real_t) const x2 = x1*x1;
    NS(be_tricub_real_t) const y2 = y1*y1;
    NS(be_tricub_real_t) const z2 = z1*z1;

    NS(be_tricub_real_t) const x3 = x2*x1;
    NS(be_tricub_real_t) const y3 = y2*y1;
    NS(be_tricub_real_t) const z3 = z2*z1;

''')

tc = ['2.0', '3.0', '4.0',  '6.0', '8.0', '9.0', '12.0', '18.0', '27.0']
def write_coef(ff,i,j,k,sp_list,jj_range,off):
    ii = i + 4 * j + 16 * k
    flag = 1
    for jj in jj_range:
        cc = tricubicMat[ii,jj]

        if cc == 0:
            continue

        if flag:
            ss = ' '
            flag = 0
        else:
            ss = '+'
                
        if abs(cc) == 1:
    #                print('coefs[%d] %s= b[%d]; '%(i,ss,j),end='')
            if ss == ' ':
                if cc == 1:
                    ff.write('    coef = b[%d];\n'%(jj-off))
                else:
                    ff.write('    coef = -b[%d];\n'%(jj-off))
            else:
                if cc == 1:
                    ss = '+'
                else:
                    ss = '-'
                ff.write('    coef %s= b[%d];\n'%(ss,jj-off))
        else:
            ind = sp_list.index(float(abs(cc)))
            if ss == ' ':
                if cc > 0:
                    ##ff.write('    coef = tri_consts[%d] * b[%d];\n'%(ind,jj))
                    ff.write('    coef = b[%d]; '%(jj-off))
                    for aa in range(abs(cc)-1):
                        ff.write('coef += b[%d]; '%(jj-off))
                    ff.write('\n')
                    #ff.write('    coef = %s * b[%d];\n'%(tc[ind],jj-off))
                else:
                    ff.write('    coef = -b[%d]; '%(jj-off))
                    for aa in range(abs(cc)-1):
                        ff.write('coef -= b[%d]; '%(jj-off))
                    ff.write('\n')
                    #ff.write('    coef = -(%s * b[%d]);\n'%(tc[ind],jj-off))
                    ##ff.write('    coef = -(tri_consts[%d] * b[%d]);\n'%(ind,jj))
            else:
                if cc > 0:
                    ss = '+'
                else:
                    ss = '-'
                #ff.write('coefs[%d] %s= b[%d]; '%(i,ss,j))
                ff.write('    ')
                for aa in range(abs(cc)-1):
                    ff.write('coef %s= b[%d]; '%(ss,jj-off))
                ff.write('\n')
                ##ff.write('    coef %s= %s * b[%d];\n'%(ss,tc[ind],jj-off))
                #ff.write('    coef %s= tri_consts[%d] * b[%d];\n'%(ss,ind,jj))
    return flag

def write_kick(ff,i,j,k):
    #ff.write('    akick = coef')
    #for _ in range(i):
    #    ff.write(' * x1')
    #for _ in range(j):
    #    ff.write(' * y1')
    #for _ in range(k):
    #    ff.write(' * z1')
    #ff.write(';\n')
    if i == 0 and j == 0 and k == 0:
        ff.write('    akick = coef;\n')
    elif i == 0 and j == 0:
        ff.write('    akick = coef * z%d;\n'%k)
    elif i == 0 and k == 0:
        ff.write('    akick = coef * y%d;\n'%j)
    elif j == 0 and k == 0:
        ff.write('    akick = coef * x%d;\n'%i)
    elif i == 0:
        ff.write('    akick = coef * (y%d * z%d);\n'%(j,k))
    elif j == 0:
        ff.write('    akick = coef * (x%d * z%d);\n'%(i,k))
    elif k == 0:
        ff.write('    akick = coef * (x%d * y%d);\n'%(i,j))
    else:
        ff.write('    akick = coef * ( (x%d * y%d) * z%d );\n'%(i,j,k))


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
#ff.write('    const double tri_consts[%d] = '%len(sp_list))
#ff.write(c_tri_list)
#ff.write(';\n\n')

def write_func(ff, sp_list, jj_range, offset):
    for i in range(4):
        for j in range(4):
            for k in range(4):
                #ff.write('    coef = coef[%d,%d,%d];\n'%(i,j,k))
                if write_coef(ff,i,j,k,sp_list,jj_range, offset): continue
                if i > 0 :
                    write_kick(ff,i-1,j,k)
                    ff.write('    ')
                    for _ in range(i):
                        ff.write('kicks[0] += akick; ')
                    ff.write('\n')
                if j > 0 :
                    write_kick(ff,i,j-1,k)
                    #ff.write('    akick = ( ( ( coef * x_power[%d] ) * y_power[%d] ) * z_power[%d] );\n'%(i,j-1,k))
                    ff.write('    ')
                    for _ in range(j):
                        ff.write('kicks[1] += akick; ')
                    ff.write('\n')
                if k > 0 :
                    write_kick(ff,i,j,k-1)
                    #ff.write('    akick = ( ( ( coef * x_power[%d] ) * y_power[%d] ) * z_power[%d] );\n'%(i,j,k-1))
                    ff.write('    ')
                    for _ in range(k):
                        ff.write('kicks[2] += akick; ')
                    ff.write('\n')

write_func(ff, sp_list, range(16), 0)
ff.write(
'''
    return;
}

SIXTRL_INLINE void NS(tricub_unscaled_kicks_from_b2)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    NS(be_tricub_real_t) x1,
    NS(be_tricub_real_t) y1,
    NS(be_tricub_real_t) z1,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* kicks
    )
{
//    NS(be_tricub_real_t) const tri_consts[9] = {2.0, 3.0, 4.0, 6.0, 8.0, 9.0, 12.0, 18.0, 27.0};
    NS(be_tricub_real_t) coef;
    NS(be_tricub_real_t) akick;

    NS(be_tricub_real_t) const x2 = x1*x1;
    NS(be_tricub_real_t) const y2 = y1*y1;
    NS(be_tricub_real_t) const z2 = z1*z1;

    NS(be_tricub_real_t) const x3 = x2*x1;
    NS(be_tricub_real_t) const y3 = y2*y1;
    NS(be_tricub_real_t) const z3 = z2*z1;

''')

write_func(ff, sp_list, range(16,32), 16)

ff.write(
'''
    return;
}

SIXTRL_INLINE void NS(tricub_unscaled_kicks_from_b3)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    NS(be_tricub_real_t) x1,
    NS(be_tricub_real_t) y1,
    NS(be_tricub_real_t) z1,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* kicks
    )
{
//    NS(be_tricub_real_t) const tri_consts[9] = {2.0, 3.0, 4.0, 6.0, 8.0, 9.0, 12.0, 18.0, 27.0};
    NS(be_tricub_real_t) coef;
    NS(be_tricub_real_t) akick;

    NS(be_tricub_real_t) const x2 = x1*x1;
    NS(be_tricub_real_t) const y2 = y1*y1;
    NS(be_tricub_real_t) const z2 = z1*z1;

    NS(be_tricub_real_t) const x3 = x2*x1;
    NS(be_tricub_real_t) const y3 = y2*y1;
    NS(be_tricub_real_t) const z3 = z2*z1;

''')

write_func(ff, sp_list, range(32,48), 32)

ff.write(
'''
    return;
}

SIXTRL_INLINE void NS(tricub_unscaled_kicks_from_b4)(
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t) const* SIXTRL_RESTRICT b,
    NS(be_tricub_real_t) x1,
    NS(be_tricub_real_t) y1,
    NS(be_tricub_real_t) z1,
    SIXTRL_ARGPTR_DEC NS(be_tricub_real_t)* kicks
    )
{
//    NS(be_tricub_real_t) const tri_consts[9] = {2.0, 3.0, 4.0, 6.0, 8.0, 9.0, 12.0, 18.0, 27.0};
    NS(be_tricub_real_t) coef;
    NS(be_tricub_real_t) akick;

    NS(be_tricub_real_t) const x2 = x1*x1;
    NS(be_tricub_real_t) const y2 = y1*y1;
    NS(be_tricub_real_t) const z2 = z1*z1;

    NS(be_tricub_real_t) const x3 = x2*x1;
    NS(be_tricub_real_t) const y3 = y2*y1;
    NS(be_tricub_real_t) const z3 = z2*z1;

''')

write_func(ff, sp_list, range(48,64), 48)

ff.write(
'''
    return;
}


#if !defined( _GPUCODE ) && defined( __cplusplus )
}
#endif /* !defined( _GPUCODE ) && defined( __cplusplus ) */

#endif /* SIXTRACKLIB_COMMON_BE_TRICUB_COEFFICIENTS_C99_H__ */
''')

