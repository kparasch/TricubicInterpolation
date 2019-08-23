#include "evaluate.h"

void tricubic_xyz_powers(double xn, double yn, double zn, double* xni, double* ynj, double* znk)
{
    xni[0] = 1.;
    ynj[0] = 1.;
    znk[0] = 1.;

    for(int i = 1; i < 4; i++)
    {
        xni[i] = xni[i-1]*xn;
        ynj[i] = ynj[i-1]*yn;
        znk[i] = znk[i-1]*zn;
    }

    return;
}

double tricubic_eval(double* coefs, double* xn, double* yn, double* zn)
{
    double result = 0.;

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                result += ( ( (coefs[i + 4 * j + 16 * k] * xn[i]) * yn[j]) * zn[k]);
            }
        }
    }

    return result;
}

double tricubic_ddx(double* coefs, double* xn, double* yn, double* zn, double dx)
{
    double result = 0.;
    

    for(int i = 1; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                result += i * ( ( (coefs[i + 4 * j + 16 * k] * xn[i-1]) * yn[j]) * zn[k]);
            }
        }
    }
    
    result /= dx;

    return result;
}


double tricubic_ddy(double* coefs, double* xn, double* yn, double* zn, double dy)
{
    double result = 0.;
    
    for(int i = 0; i < 4; i++)
    {
        for(int j = 1; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                result += j * ( ( (coefs[i + 4 * j + 16 * k] * xn[i]) * yn[j-1]) * zn[k]);
            }
        }
    }

    result /= dy;

    return result;
}

double tricubic_ddz(double* coefs, double* xn, double* yn, double* zn, double dz)
{
    double result = 0.;

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 1; k < 4; k++)
            {
                result += k * ( ( (coefs[i + 4 * j + 16 * k] * xn[i]) * yn[j]) * zn[k-1]);
            }
        }
    }

    result /= dz;

    return result;
}


double tricubic_ddxdy(double* coefs, double* xn, double* yn, double* zn, double dx, double dy)
{
    double result = 0.;
    
    for(int i = 1; i < 4; i++)
    {
        for(int j = 1; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                result += (i * j) * ( ( (coefs[i + 4 * j + 16 * k] * xn[i-1]) * yn[j-1]) * zn[k]);
            }
        }
    }

    result /= (dx * dy);

    return result;
}

double tricubic_ddxdz(double* coefs, double* xn, double* yn, double* zn, double dx, double dz)
{
    double result = 0.;
    
    for(int i = 1; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 1; k < 4; k++)
            {
                result += (i * k) * ( ( (coefs[i + 4 * j + 16 * k] * xn[i-1]) * yn[j]) * zn[k-1]);
            }
        }
    }

    result /= (dx * dz);

    return result;
}

double tricubic_ddydz(double* coefs, double* xn, double* yn, double* zn, double dy, double dz)
{
    double result = 0.;
    
    for(int i = 0; i < 4; i++)
    {
        for(int j = 1; j < 4; j++)
        {
            for(int k = 1; k < 4; k++)
            {
                result += (j * k) * ( ( (coefs[i + 4 * j + 16 * k] * xn[i]) * yn[j-1]) * zn[k]);
            }
        }
    }

    result /= (dy * dz);

    return result;
}

double tricubic_ddxdydz(double* coefs, double* xn, double* yn, double* zn, double dx, double dy, double dz)
{
    double result = 0.;
    
    for(int i = 1; i < 4; i++)
    {
        for(int j = 1; j < 4; j++)
        {
            for(int k = 1; k < 4; k++)
            {
                result += ( (i * j) * k) * ( ( (coefs[i + 4 * j + 16 * k] * xn[i-1]) * yn[j-1]) * zn[k-1]);
            }
        }
    }

    result /= ((dx * dy) * dz);

    return result;
}

double tricubic_ddx2(double* coefs, double* xn, double* yn, double* zn, double dx)
{
    double result = 0.;
    

    for(int i = 2; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                result += ( i * (i - 1) ) * ( ( (coefs[i + 4 * j + 16 * k] * xn[i-2]) * yn[j]) * zn[k]);
            }
        }
    }
    
    result /= dx*dx;

    return result;
}

double tricubic_ddy2(double* coefs, double* xn, double* yn, double* zn, double dy)
{
    double result = 0.;
    
    for(int i = 0; i < 4; i++)
    {
        for(int j = 2; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                result += ( j * (j - 1) ) * ( ( (coefs[i + 4 * j + 16 * k] * xn[i]) * yn[j-2]) * zn[k]);
            }
        }
    }

    result /= dy*dy;

    return result;
}

double tricubic_ddz2(double* coefs, double* xn, double* yn, double* zn, double dz)
{
    double result = 0.;

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 2; k < 4; k++)
            {
                result += ( k * (k - 1) ) * ( ( (coefs[i + 4 * j + 16 * k] * xn[i]) * yn[j]) * zn[k-2]);
            }
        }
    }

    result /= dz*dz;

    return result;
}

