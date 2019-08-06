#include "evaluate.h"

void xyz_powers(double xn, double yn, double zn, double* xni, double* ynj, double* znk)
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

double eval(double* coefs, double* xn, double* yn, double* zn)
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

double ddx(double* coefs, double* xn, double* yn, double* zn, double dx)
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


double ddy(double* coefs, double* xn, double* yn, double* zn, double dy)
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

double ddz(double* coefs, double* xn, double* yn, double* zn, double dz)
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


double ddxdy(double* coefs, double* xn, double* yn, double* zn, double dx, double dy)
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

double ddxdz(double* coefs, double* xn, double* yn, double* zn, double dx, double dz)
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

double ddydz(double* coefs, double* xn, double* yn, double* zn, double dy, double dz)
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

double ddxdydz(double* coefs, double* xn, double* yn, double* zn, double dx, double dy, double dz)
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
