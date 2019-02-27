/*
 * nbody_psm_serial.c
 *
 *
 * Serialized port of the NBody PSM method
 * TODO: Documentation
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include "nbody_psm.h"

/*
 * Calculates the Cauchy Procuct of a vector 'a'  raised to a 
 * power 'pow' and stored it in vector 'b'.
 */
double cauchy_power(double a[], double b[], int n, double pow) 
{
    double cauchy_pow = 0;
    for (int i = 0; i < n - 1; i++)
    {
        cauchy_pow += (i * (pow + 1) + pow - n) * a[i + 1] * b[n - i];
    }
    return cauchy_pow;
}

/*
 * Calculates the Cauchy Procuct of two vectors a and b.
 */
double cauchy_prod(double a[], double b[], int n)
{
    double prod = 0;

    for (int i = 0; i < n; i++)
    {
        prod += a[i] * b[n - i];
    }

    return prod;
}

/*
 * Calculates the Cauchy Product of a vector 'a' squared.
 */
double cauchy_square(double a[], int n)
{
    double square = 0;

    int limit = (n/2) - 1;

    if (n % 2 == 0) 
    {
        for (int i = 0; i <= limit; i++) 
        {
            square += a[i] * a[n - i];
        }
        square = (2 * square) + (a[n/2] * a[n/2]);
    }
    return square;
}

/*
 * Calculates the value of a polynomial at a given 't' 
 * using Horner's algorithm.
 */
double horner_value(double c[], double t, int n)
{
    int sum = c[n - 1] + (t * c[n]);
    for (int i = n - 2; i >= 0; i--)
    {
        sum = c[i] + t * sum;
    }
    return sum;
}

 int main() 
 {
    return 0;
 }
