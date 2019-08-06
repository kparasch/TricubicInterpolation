#include "coords.h"

void tricubic_coords_to_indices(double x, double y, double z, double x0, double y0, double z0, double dx, double dy, double dz, int &ix, int &iy, int &iz)
{
    double fx = (x - x0)/dx;
    double fy = (y - y0)/dy;
    double fz = (z - z0)/dz;

    ix = (int)fx;
    iy = (int)fy;
    iz = (int)fz;

    return;
}

int tricubic_coords_to_indices_and_floats(double x, double y, double z, double x0, double y0, double z0, double dx, double dy, double dz, int &ix, int &iy, int &iz, double &xn, double &yn, double &zn, int ix_bound_low, int ix_bound_up, int iy_bound_low, int iy_bound_up, int iz_bound_low, int iz_bound_up)
{


    double fx = (x - x0)/dx;
    double fy = (y - y0)/dy;
    double fz = (z - z0)/dz;

    ix = (int)fx;
    iy = (int)fy;
    iz = (int)fz;

    xn = fx - ix;
    yn = fy - iy;
    zn = fz - iz;

    int inside_box = 1;
    if( ix < ix_bound_low || ix > ix_bound_up )
        inside_box = 0;
    else if( iy < iy_bound_low || iy > iy_bound_up )
        inside_box = 0;
    else if( iz < iz_bound_low || iz > iz_bound_up )
        inside_box = 0;

    return inside_box;

}
