#ifndef __DERIVS_H__
#define __DERIVS_H__

double* finite_diff(int shape1, int shape2, double A[][shape1][shape2], int ix, int iy, int iz);
double* exact_diff(int shape1, int shape2, int shape3, double A[][shape1][shape2][shape3], int ix, int iy, int iz, double dx, double dy, double dz);

#endif
