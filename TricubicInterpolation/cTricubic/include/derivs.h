#ifndef __DERIVS_H__
#define __DERIVS_H__

#include <stdlib.h>

double* tricubic_finite_diff(int shape1, int shape2, double* A, int ix, int iy, int iz);
double* tricubic_exact_diff(int shape1, int shape2, int shape3, double* A, int ix, int iy, int iz, double dx, double dy, double dz);

#endif
