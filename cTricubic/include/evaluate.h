#ifndef __EVALUATE_H__
#define __EVALUATE_H__

double xyz_powers(double xn, double yn, double zn, double* xni, double* ynj, double* znk);
double eval(double* coefs, double* xn, double* yn, double* zn);
double ddx(double* coefs, double* xn, double* yn, double* zn, double dx);
double ddy(double* coefs, double* xn, double* yn, double* zn, double dy);
double ddz(double* coefs, double* xn, double* yn, double* zn, double dz);
double ddxdy(double* coefs, double* xn, double* yn, double* zn, double dx, double dy);
double ddxdz(double* coefs, double* xn, double* yn, double* zn, double dx, double dz);
double ddydz(double* coefs, double* xn, double* yn, double* zn, double dy, double dz);
double ddxdydz(double* coefs, double* xn, double* yn, double* zn, double dx, double dy, double dz);

#endif
