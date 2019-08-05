#include "evaluate.h"

double eval(double* coefs, double xn, double yn, double zn)
{
    double result = 0.;
    
    double xni = 1./xn; //xni = xn**i
    double ynj = 1./yn; //ynj = yn**j
    double znk = 1./zn; //znk = zn**k

    for(int i = 0; i < 4; i++)
    {
        xni *= xn; 
        for(int j = 0; j < 4; j++)
        {
            ynj *= yn;
            for(int k = 0; k < 4; k++)
            {
                znk *= zn;
                result += ( ( (coefs[i + 4 * j + 16 * k] * xni) * ynj) * znk);
            }
        }
    }

    return result;
}

double ddx(double* coefs, double xn, double yn, double zn, double dx)
{
    double result = 0.;
    
    double xni = 1./xn; //xni = xn**i
    double ynj = 1./yn; //ynj = yn**j
    double znk = 1./zn; //znk = zn**k

    for(int i = 1; i < 4; i++)
    {
        xni *= xn; 
        for(int j = 0; j < 4; j++)
        {
            ynj *= yn;
            for(int k = 0; k < 4; k++)
            {
                znk *= zn;
                result += i * ( ( (coefs[i + 4 * j + 16 * k] * xni) * ynj) * znk);
            }
        }
    }
    
    result /= dx;

    return result;
}


double ddy(double* coefs, double xn, double yn, double zn, double dy)
{
    double result = 0.;
    
    double xni = 1./xn; //xni = xn**i
    double ynj = 1./yn; //ynj = yn**j
    double znk = 1./zn; //znk = zn**k

    for(int i = 0; i < 4; i++)
    {
        xni *= xn; 
        for(int j = 1; j < 4; j++)
        {
            ynj *= yn;
            for(int k = 0; k < 4; k++)
            {
                znk *= zn;
                result += j * ( ( (coefs[i + 4 * j + 16 * k] * xni) * ynj) * znk);
            }
        }
    }

    result /= dy;

    return result;
}

double ddz(double* coefs, double xn, double yn, double zn, double dz)
{
    double result = 0.;
    
    double xni = 1./xn; //xni = xn**i
    double ynj = 1./yn; //ynj = yn**j
    double znk = 1./zn; //znk = zn**k

    for(int i = 0; i < 4; i++)
    {
        xni *= xn; 
        for(int j = 0; j < 4; j++)
        {
            ynj *= yn;
            for(int k = 1; k < 4; k++)
            {
                znk *= zn;
                result += k * ( ( (coefs[i + 4 * j + 16 * k] * xni) * ynj) * znk);
            }
        }
    }

    result /= dz;

    return result;
}


double ddxdy(double* coefs, double xn, double yn, double zn, double dx, double dy)
{
    double result = 0.;
    
    double xni = 1./xn; //xni = xn**i
    double ynj = 1./yn; //ynj = yn**j
    double znk = 1./zn; //znk = zn**k

    for(int i = 1; i < 4; i++)
    {
        xni *= xn; 
        for(int j = 1; j < 4; j++)
        {
            ynj *= yn;
            for(int k = 0; k < 4; k++)
            {
                znk *= zn;
                result += (i * j) * ( ( (coefs[i + 4 * j + 16 * k] * xni) * ynj) * znk);
            }
        }
    }

    result /= (dx * dy);

    return result;
}

double ddxdz(double* coefs, double xn, double yn, double zn, double dx, double dz)
{
    double result = 0.;
    
    double xni = 1./xn; //xni = xn**i
    double ynj = 1./yn; //ynj = yn**j
    double znk = 1./zn; //znk = zn**k

    for(int i = 1; i < 4; i++)
    {
        xni *= xn; 
        for(int j = 0; j < 4; j++)
        {
            ynj *= yn;
            for(int k = 1; k < 4; k++)
            {
                znk *= zn;
                result += (i * k) * ( ( (coefs[i + 4 * j + 16 * k] * xni) * ynj) * znk);
            }
        }
    }

    result /= (dx * dz);

    return result;
}

double ddydz(double* coefs, double xn, double yn, double zn, double dy, double dz)
{
    double result = 0.;
    
    double xni = 1./xn; //xni = xn**i
    double ynj = 1./yn; //ynj = yn**j
    double znk = 1./zn; //znk = zn**k

    for(int i = 0; i < 4; i++)
    {
        xni *= xn; 
        for(int j = 1; j < 4; j++)
        {
            ynj *= yn;
            for(int k = 1; k < 4; k++)
            {
                znk *= zn;
                result += (j * k) * ( ( (coefs[i + 4 * j + 16 * k] * xni) * ynj) * znk);
            }
        }
    }

    result /= (dy * dz);

    return result;
}

double ddxdydz(double* coefs, double xn, double yn, double zn, double dx, double dy, double dz)
{
    double result = 0.;
    
    double xni = 1./xn; //xni = xn**i
    double ynj = 1./yn; //ynj = yn**j
    double znk = 1./zn; //znk = zn**k

    for(int i = 1; i < 4; i++)
    {
        xni *= xn; 
        for(int j = 1; j < 4; j++)
        {
            ynj *= yn;
            for(int k = 1; k < 4; k++)
            {
                znk *= zn;
                result += ( (i * j) * k) * ( ( (coefs[i + 4 * j + 16 * k] * xni) * ynj) * znk);
            }
        }
    }

    result /= ((dx * dy) * dz);

    return result;
}
