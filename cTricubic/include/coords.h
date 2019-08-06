#ifndef __TRICUBIC_COORDS_H__
#define __TRICUBIC_COORDS_H__

void tricubic_coords_to_indices(double x, double y, double z, double x0, double y0, double z0, double dx, double dy, double dz, int &ix, int &iy, int &iz);
int tricubic_coords_to_indices_and_floats(double x, double y, double z, double x0, double y0, double z0, double dx, double dy, double dz, int &ix, int &iy, int &iz, double &xn, double &yn, double &zn, int ix_bound_low, int ix_bound_up, int iy_bound_low, int iy_bound_up, int iz_bound_low, int iz_bound_up);

#endif
