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
#include <unistd.h>
#include <getopt.h>
#include "nbody_psm.h"

int time_step;                  // Time step
int num_ts;                     // Number of time steps
int granularity;                // Output granularity (in time steps)
int mac_degree = MAC_DEGREE;    // Degree of maclaurin polynomials
int num_bodies = NUM_BODIES;    // Number of bodies

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

/*
 * Prints the usage
 */
void usage(int argc, char *argv[])
{
    printf("Usage: %s [-tng] <input_file>\n", argv[0]);
    printf("\t-t: Specify the time step")
    printf("\t-n: Specify the number of time steps to complete")
    printf("\t-g: Specify the granularity of the output (how often to report new state)")
}

int main(int argc, char *argv[]) 
{
    // Get command line arguments

    int c;
    while ((c = getopt(argc, argv, "t:n:g:")) != -1) 
    {
        switch (c) 
        {
        case 'n':
            num_ts = atoi(optarg);
            break;
        case 'g':
            granularity = atoi(optarg);
            break;
        case 't':
            time_step = atoi(optarg);
            break;
        default:
            usage(argc, argv);
            exit(EXIT_FAILURE);
        }
    }
    if (optind != argc - 1) {
        usage(argc, argv);
        exit(EXIT_FAILURE);
    }

    // Read input file
    

    // Body positions (in sequence)
    double x[num_bodies + 1][mac_degree + 1], y[num_bodies + 1][mac_degree + 1], y[num_bodies + 1][mac_degree + 1];
    // Body velocities (in sequence)
    double u[num_bodies + 1][mac_degree + 1], v[num_bodies + 1][mac_degree + 1], w[num_bodies + 1][mac_degree + 1];
    // Body masses
    double mass[num_bodies + 1]
    // Used for calculating body values
    double r[num_bodies + 1][num_bodies + 1][mac_degree + 1], b[num_bodies + 1][num_bodies + 1][mac_degree + 1]
    // PSM values for updating body position and velocity
    double x_psm, y_psm, z_psm, u_psm, v_psm, w_psm;



    // Cleanup and exit
    free(fname)

    return 0;
}
