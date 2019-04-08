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
#include <omp.h>

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
#   ifdef _OPENMP
#   pragma omp parallel for
#   endif
    for (int i = n - 2; i >= 0; i--)
    {
        sum = c[i] + t * sum;
    }
    return sum;
}

double perform_calculation(FILE *in_file)
{
    struct timeval tv;

    if (debug)
    {
        printf("Initializing Body Members\n");
    }
    
    if (debug)
    {
        printf("mac_degree: %d\n", mac_degree);
    }
    

    // Body positions (in sequence)
    double x[num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("x Initialized\n");
    }
    double y[num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("y Initialized\n");
    }
    double z[num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("z Initialized\n");
    }
    
    if (debug)
    {
        printf("Positions Initialized\n");
    }
    
    // Body velocities (in sequence)
    double u[num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("u Initialized\n");
    }
    double v[num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("v Initialized\n");
    }
    double w[num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("w Initialized\n");
    }
    
    // Body masses
    double mass[num_bodies + 1];
    if (debug)
    {
        printf("mass Initialized\n");
    }
    
    // Used for calculating body values
    double X[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("X Initialized\n");
    }
    double Y[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("Y Initialized\n");
    }
    double Z[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("Z Initialized\n");
    }
    double r[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("r Initialized\n");
    }
    double b[num_bodies + 1][num_bodies + 1][mac_degree + 1];
    if (debug)
    {
        printf("b Initialized\n");
    }
    
    // PSM values for updating body position and velocity
    double x_psm, y_psm, z_psm, u_psm, v_psm, w_psm;
    
    if (debug)
    {
        printf("Body Members Initialized\n");
    }
    
    // File pointers for body data output
    FILE *body_output[num_bodies + 1];
    char outfile_name[16];
    char line[BUFFER_SIZE];
    
    if (debug)
    {
        printf("Starting file import\n");
    }
    
    double tmp[7];
    for (int i = 1; i <= num_bodies; i++)
    {
        if (debug)
        {
            printf("Reading body %d\n", i);
        }

        for (int j = 0; j < 7; j++)
        {
            if (fgets(line, BUFFER_SIZE, in_file) == NULL)
            {
                printf("Error reading input file\n");
                exit(1);
            }
            tmp[j] = atof(line);  

        }
        
        mass[i] = tmp[0];
        x[i][0] = tmp[1];
        y[i][0] = tmp[2];
        z[i][0] = tmp[3];
        u[i][0] = tmp[4];
        v[i][0] = tmp[5];
        w[i][0] = tmp[6];

        printf("Body %d\n", i);
        printf("Mass: %f\n", mass[i]);
        printf("x: %f\ty: %f\tz: %f\n", x[i][0], y[i][0], z[i][0]);
        printf("u: %f\tv: %f\tw: %f\n", u[i][0], v[i][0], w[i][0]);
        printf("\n");
        
    }

    // TODO: Replace hard coded values with values read from input file

	/*    Set up the initial positions and velocities  */
    if (verbose) 
    {
        printf("Initial Body States:\n");   
        for (int i = 1; i <= 10; i++)
        {
            printf("Body %d\n", i);
            printf("Mass: %f\n", mass[i]);
            printf("x: %f\ty: %f\tz: %f\n", x[i][0], y[i][0], z[i][0]);
            printf("u: %f\tv: %f\tw: %f\n", u[i][0], v[i][0], w[i][0]);
            printf("\n");
        }
    }
    
    if (output)
    {
        if (debug)
        {
            printf("Opening files for output\n\n");
        }
        
        for (int i = 0; i < num_bodies; i++)
        {
            sprintf(outfile_name, "body%04d", i+1);
            if ((body_output[i + 1] = fopen(outfile_name, "w")) == NULL)
            {
                printf("Output File %s could not be opened for writing: %s\n", outfile_name,
                    strerror(errno));
            }
        }
    }
    
    if (verbose)
    {
        printf("Time step: %f\n", time_step);
        printf("Number of time steps: %d\n", num_ts);
        printf("Granularity: %d\n\n", granularity);
    }

    if (debug)
    {
        printf("Starting Main Loop with %d time steps\n", num_ts);
    }
    
    int step, i, j, k;

    START_TIMER(psm)

    for (step = 0; step <= num_ts; step++)
    {
        if (debug || verbose)
        {
            if (step % granularity == 0)
                printf("Step %d\n", step);
        }

        // Begin parallel region
#       ifdef _OPENMP 
#       pragma omp parallel
#       endif
        {   
            // LOOP 1
#           ifdef _OPENMP
#           pragma omp for
#           endif
            for (i = 1; i <= num_bodies; i++)
            {
                // LOOP 2
                for(j = 1; j <= i - 1; j++)
                {
                    X[i][j][0] = x[j][0] - x[i][0];
                    Y[i][j][0] = y[j][0] - y[i][0];
                    Z[i][j][0] = z[j][0] - z[i][0];
                    r[i][j][0] = pow(X[i][j][0], 2) + pow(Y[i][j][0], 2) + pow(Z[i][j][0], 2);
                    b[i][j][0] = pow(r[i][j][0], -1.5);
                }

                // LOOP 3
                for (j = i + 1; j < num_bodies; j++)
                {
                    X[i][j][0] = x[j][0] - x[i][0];
                    Y[i][j][0] = y[j][0] - y[i][0];
                    Z[i][j][0] = z[j][0] - z[i][0];
                    r[i][j][0] = pow(X[i][j][0], 2) + pow(Y[i][j][0], 2) + pow(Z[i][j][0], 2);
                    b[i][j][0] = pow(r[i][j][0], -1.5);
                }
            }
            
#           ifdef _OPENMP
#           pragma omp barrier
#           endif

            // LOOP 4
#           ifdef _OPENMP
#           pragma omp for
#           endif
            for (k = 1; k <= mac_degree; k++)
            {
                // LOOP 5              
                for (i = 1; i <= num_bodies; i++)
                {
                    x[i][k] = u[i][k - 1] / k;
                    y[i][k] = v[i][k - 1] / k;
                    z[i][k] = w[i][k - 1] / k;
                }

                // LOOP 6
                for (i = 1; i <= num_bodies; i++)
                {
                    // LOOP 7
#                   ifdef _OPENMP
#                   pragma omp single
#                   endif                   
                    {
                        for (j = 1; j <= num_bodies; j++)
                        {
                            X[i][j][k] = x[j][k] - x[i][k];
                            Y[i][j][k] = y[j][k] - y[i][k];
                            Z[i][j][k] = z[j][k] - z[i][k];
                            r[i][j][k] = cauchy_prod(X[i][j], X[i][j], k) + 
                                            cauchy_prod(Y[i][j], Y[i][j], k) + 
                                            cauchy_prod(Z[i][j], Z[i][j], k);
                            b[i][j][k] = cauchy_power(r[i][j], b[i][j], k - 1, -1.5);
                        }
                    }
                    

                    // LOOP 8
#                   ifdef _OPENMP
#                   pragma omp single
#                   endif                   
                    {
                        for (j = i + 1; j <= num_bodies; j++)
                        {
                            X[i][j][k] = x[j][k] - x[i][k];
                            Y[i][j][k] = y[j][k] - y[i][k];
                            Z[i][j][k] = z[j][k] - z[i][k];
                            r[i][j][k] = cauchy_prod(X[i][j], X[i][j], k) + 
                                            cauchy_prod(Y[i][j], Y[i][j], k) + 
                                            cauchy_prod(Z[i][j], Z[i][j], k);
                            b[i][j][k] = cauchy_power(r[i][j], b[i][j], k - 1, -1.5);
                        }
                    }
                }
                
#               ifdef _OPENMP
#               pragma omp barrier
#               endif   
                
                // LOOP 9
                for (i = 1; i <= num_bodies; i++)
                {
                    u[i][k] = 0;
                    v[i][k] = 0;
                    w[i][k] = 0;
                    
                    // LOOP 10
                    for (j = 1; j <= i - 1; j++)
                    {
                        u[i][k] += mass[j] * cauchy_prod(X[i][j], b[i][j], k - 1);
                        v[i][k] += mass[j] * cauchy_prod(Y[i][j], b[i][j], k - 1);
                        w[i][k] += mass[j] * cauchy_prod(Z[i][j], b[i][j], k - 1);
                    }
                    
                    // LOOP 11
                    for (int j = i + 1; j <= num_bodies; j++)
                    {
                        u[i][k] += mass[j] * cauchy_prod(X[i][j], b[i][j], k - 1);
                        v[i][k] += mass[j] * cauchy_prod(Y[i][j], b[i][j], k - 1);
                        w[i][k] += mass[j] * cauchy_prod(Z[i][j], b[i][j], k - 1);
                    }
                    
                    u[i][k] = u[i][k] / k;
                    v[i][k] = v[i][k] / k;
                    w[i][k] = w[i][k] / k;
                }
            }
        }

        // Determine the values of the Maclaurin polynomial using Horner's algorithm and the
        // stored Maclauren coefficients
                        
#       ifdef _OPENMP
#       pragma omp barrier
#       pragma omp parallel
#       endif   
        
        {
#           ifdef _OPENMP
#           pragma omp for
#           endif 
            for (i = 1; i <= num_bodies; i++)
            {
                x_psm = horner_value(x[i], time_step, mac_degree);
                x[i][0] = x_psm;
                
                y_psm = horner_value(y[i], time_step, mac_degree);
                y[i][0] = y_psm;
                
                z_psm = horner_value(z[i], time_step, mac_degree);
                z[i][0] = z_psm;
                
                u_psm = horner_value(u[i], time_step, mac_degree);
                u[i][0] = u_psm;
                
                v_psm = horner_value(v[i], time_step, mac_degree);
                v[i][0] = v_psm;
                
                w_psm = horner_value(w[i], time_step, mac_degree);
                w[i][0] = w_psm;
            }
        }
        

        // Output the step number based on the output granularity
        if (step % granularity == 0)
        {   
            if (debug) 
            {
                for (i = 1; i <= num_bodies; i++)
                {
                    printf("body %d: x: %lf\ty: %lf\tz: %lf\n", i, x[i][0], y[i][0], z[i][0]);
                }
            }
        }
        
        if (output)
        {
            for (i = 0; i < num_bodies; i++)
            {
                fprintf(body_output[i + 1], "x: %lf\ty: %lf\tz: %lf\n", x[i][0], y[i][0], z[i][0]);
            }
        }
    }
    STOP_TIMER(psm)
    return GET_TIMER(psm);

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
    double time;

    
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
    
    if (debug)
    {
        printf("Starting calculation\n");
    }
    time = perform_calculation(in_file);
    
    if (timer)
    {
        printf("Bodies: %d\tCalculation: %0.4fs\n", num_bodies, time);
    }
    
    if (debug)
    {
        printf("End of calculation\n\n");
    }
    
    return 0;
}
