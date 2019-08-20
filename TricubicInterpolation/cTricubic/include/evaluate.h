#ifndef __EVALUATE_H__
#define __EVALUATE_H__

void tricubic_xyz_powers(double xn, double yn, double zn, double* xni, double* ynj, double* znk);
double tricubic_eval(double* coefs, double* xn, double* yn, double* zn);
double tricubic_ddx(double* coefs, double* xn, double* yn, double* zn, double dx);
double tricubic_ddy(double* coefs, double* xn, double* yn, double* zn, double dy);
double tricubic_ddz(double* coefs, double* xn, double* yn, double* zn, double dz);
double tricubic_ddxdy(double* coefs, double* xn, double* yn, double* zn, double dx, double dy);
double tricubic_ddxdz(double* coefs, double* xn, double* yn, double* zn, double dx, double dz);
double tricubic_ddydz(double* coefs, double* xn, double* yn, double* zn, double dy, double dz);
double tricubic_ddxdydz(double* coefs, double* xn, double* yn, double* zn, double dx, double dy, double dz);

#endif
