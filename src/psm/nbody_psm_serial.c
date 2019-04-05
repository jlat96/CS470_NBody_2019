/*
 * nbody_psm_serial.c
 *
 *
 * Serialized port of the NBody PSM method
 *
 * Input:
 *
 * An input file begins with an integer value denoting the number of 
 * bodies in the system The next 7 values are the body's mass, the x, 
 * y, and z values for its position relative to the center of the system,
 * and the body's u, v , and w velocity values.
 *
 * A sample file with 2 bodies would appear as follows
 *
 * 2
 * 1
 * 0
 * 0
 * 0
 * 0
 * 0
 * 0
 * 7.6923076923076926E-9
 * -12.10226300000
 * -26.73256000000
 * 6.362842900000
 * 0.1714028159354
 * -.1021868903979
 * -.3854629379438E-01
 *
 *
 * Acknowledgements:
 *
 * Timing macros provided by Dr. Mike Lam
 * PSM Algorithm based on work and code written by Dr. Sochacki
 * CS470 Project from 2017 about the n-body problem was viewed as
 * reference when porting Dr. Sochacki's Matlab code (Provided by Dr. Lam)
 *
 * Roadmap for next steps:
 * 1) Ensure that results from the serial version are accurate (Week 2/3)
 * 2) Change variables used in calculation to exist at file scope (dynamically allocated memory) (Week 3)
 * 3) Decompose steps in PSM algorithm to functions (Week 3/4)
 * 4) Assess parallelizability of functions/regions with OpenMP (Week 4)
 * 5) OpenMP Analysis (Week 4/5)
 * 6) Analysis comparison with other N-Body algorithms (Week 5)
 * 
 * TODO: Analysis and detailed documentation
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/time.h>

#include "nbody_psm.h"

#define BUFFER_SIZE 22

double time_step = DEFAULT_TS;          // Time step
int num_ts = DEFAULT_NUM;               // Number of time steps
int granularity = DEFAULT_GRANULARITY;  // Output granularity (in time steps)
int mac_degree = MAC_DEGREE;            // Degree of maclaurin polynomials
int num_bodies = NUM_BODIES;            // Number of bodies

bool debug = false;
bool verbose = false;
bool timer = false;
bool output = false;

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
 * Calculates the value of a polynomial at a given 't' 
 * using Horner's algorithm.
 */
double horner_value(double c[], double t, int n)
{
    double sum = c[n - 1] + (t * c[n]);
    for (int i = n - 2; i >= 0; i--)
    {
        sum = c[i] + t * sum;
    }
    return sum;
}

void load_bodies(int num_bodies, FILE *in_file, double *mass, double *x, 
                    double *y, double *z, double *u, double *v, double *w)
{
    char line[BUFFER_SIZE];
    double tmp[7];
    for (int i = 1; i <= num_bodies; i++)
    {
        printf("Reading body %d\n", i);
        for (int j = 0; j < 7; j++)
        {
            if (fgets(line, BUFFER_SIZE, in_file) == NULL)
            {
                printf("Error reading input file\n");
                return 1;
            }
            tmp[j] = atof(line);  
            printf("%f\n", tmp[j]);      
        }


        
        
        
    }
    
    if (debug)
    {
        printf("File read\n");
    }
    
    printf("Step %d\n", step);
    for (int i = 1; i <= num_bodies; i++)
    {
        //printf("body %d: x: %lf\ty: %lf\tz: %lf\n", i, x[i][0], y[i][0], z[i][0]);
    }
}

void perform_calculation(FILE *in_file)
{
// Body positions (in sequence)
    double x[num_bodies + 1][mac_degree + 1];
    double y[num_bodies + 1][mac_degree + 1];
    double z[num_bodies + 1][mac_degree + 1];
    
    // Body velocities (in sequence)
    double u[num_bodies + 1][mac_degree + 1];
    double v[num_bodies + 1][mac_degree + 1];
    double w[num_bodies + 1][mac_degree + 1];
    
    // Body masses
    double mass[num_bodies + 1];
    
    // Used for calculating body values
    double X[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    double Y[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    double Z[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    double r[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    double b[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    
    // PSM values for updating body position and velocity
    double x_psm, y_psm, z_psm, u_psm, v_psm, w_psm;
    
    // File pointers for body data output
    FILE *body_output[num_bodies + 1];
    char outfile_name[16];
    
    load_bodies(num_bodies, in_file, mass, x, y, z, u, v, w);
    
    if (debug)
    {
        printf("Bodies Loaded\n");
    }
}

/*
 * Prints the usage
 */
void usage(char *argv[])
{
    printf("Usage: %s [-dghntTv]\n", argv[0]);
    printf("\t-d: Print debug output\n");
    printf("\t-g: Specify the granularity of the output (how often to report new state)\n");
    printf("\t-h: Print usage\n");
    printf("\t-n: Specify the number of time steps to complete\n");
    printf("\t-t: Specify the time step\n");
    printf("\t-T: Print timer output\n");
    printf("\t-v: Print verbose output\n");

}

int main(int argc, char *argv[]) 
{
    // Get command line arguments
    int c;
    
    struct timeval tv;
    
    while ((c = getopt(argc, argv, "dg:hn:ot:Tv")) != -1) 
    {
        switch (c) 
        {
        case 'd':
            debug = true;
            break;
        case 'g':
            granularity = atoi(optarg);
            break;
        case 'h':
            usage(argv);
            exit(EXIT_FAILURE);    
        case 'o':
            output = true;
            break;   
        case 'n':
            num_ts = atoi(optarg);
            break;
        case 't':
            time_step = atoi(optarg);
            break;
        case 'T':
            timer = true;
            break;
        case 'v':
            verbose = true;
            break;
        default:
            usage(argv);
            exit(EXIT_FAILURE);
        }
    }

    if (optind != argc - 1)
    {
        usage(argv);
        exit(EXIT_FAILURE);
    } 
    
    FILE *in_file;
    char line[BUFFER_SIZE];
    char *file_name = argv[optind];

    if (debug)
    {
        printf("filename: %s\n", file_name);
    }
    
    in_file = fopen(file_name, "r");
    if (in_file == NULL) 
    {
        printf("The file %s could not be opened!\n", file_name);
        return 1;
    }

    // Read input file, determine number of bodies, initialize arrays to size
    if (fgets(line, BUFFER_SIZE, in_file) == NULL)
    {
        printf("Error reading the input file\n");
        return 1;
    }

    num_bodies = atoi(line);
    
    if (debug)
    {
        printf("Num Bodies: %d\n", num_bodies);  
    }
    
    perform_calculation(in_file);
    
    return 0;
}
