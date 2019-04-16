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

/*
 * Calculates the Cauchy Procuct of a vector 'a'  raised to a 
 * power 'pow' and stored it in vector 'b'.
 */
double cauchy_power(double a[], double b[], int n, double pow) 
{
    double cauchy_pow = 0;
    //Not Parallelizable
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
    
    //Not Parallelizable
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
	
    //Parallelizable
    for (int i = n - 2; i >= 0; i--)
    {
        sum = c[i] + t * sum;
    }
    return sum;
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
    bool debug = false;
    bool verbose = false;
    bool timer = false;
    
    struct timeval tv;
    
    while ((c = getopt(argc, argv, "dg:hn:t:Tv")) != -1) 
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
    
    /*
    if (optind != argc - 1) {
        usage(argv);
        exit(EXIT_FAILURE);
    } 
    */   


    // TODO: Fix File Input

    /*
    FILE *in_file;
    char line[BUFFER_SIZE];
    char *file_name = argv[optind];

    printf("filename: %s\n", file_name);

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
    */
    
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

    /*
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
            tmp[i] = atof(line);  
            printf("%f\n", tmp[i]);      
        }

        for (int j = 0; j < 7; j++)
        {
            printf("%f\n", tmp[i]);
        }

        mass[i] = tmp[0];
        x[i][0] = tmp[1];
        y[i][0] = tmp[2];
        z[i][0] = tmp[3];
        u[i][0] = tmp[4];
        v[i][0] = tmp[5];
        w[i][0] = tmp[6];
    }
    */
    

    // TODO: Replace hard coded valyes with values read from input file
    
    mass[1] =  1.;
    mass[2] =  1.66013679527193035E-7;
    mass[3] =  2.44783959796682464E-6;
    mass[4] =  3.04043273871083524E-6;
    mass[5] =  3.22714936215392876E-7;
    mass[6] =  9.54790662147324233E-4;
    mass[7] =  2.85877644368210402E-4;
    mass[8] =  4.35540069686411149E-5;
    mass[9] =  5.17759138448793649E-5;
    mass[10] = 7.6923076923076926E-9;

	/*    Set up the initial positions and velocities  */
    /* SUN */
	
    //Now to be done via input file
    x[1][0] = 0;
  	y[1][0] = 0;
  	z[1][0] = 0;
  	u[1][0] = 0;
  	v[1][0] = 0;
  	w[1][0] = 0;

    /* MERCURY */
  	x[2][0] = 0.1633519000000E-02;
  	y[2][0] = 0.3077199200000;
  	z[2][0] = 0.2498653800000E-01;
  	u[2][0] = -1.963534107001;
  	v[2][0] = 0.6841909835729E-01;
  	w[2][0] = 0.1858224458145;

    /* VENUS */
	x[3][0] = 0.2640622900000;
  	y[3][0] = 0.6709840300000;
  	z[3][0] = -.6078524800000E-02;
  	u[3][0] = -1.097989091627;
  	v[3][0] = 0.4250727728824;
  	w[3][0] = 0.6918653377540E-01;
	
    /* EARTH */
	x[4][0] =  0.6643154100000E-01;
  	y[4][0] =  0.9817401900000;
  	z[4][0] =  0.6625301000000E-05;
  	u[4][0] =  -1.014141067952;
  	v[4][0] =  0.6377084582524E-01;
  	w[4][0] =  -.1047343992878E-05;

    /* MARS */ 
	x[5][0] = 1.105892500000;
  	y[5][0] = -.8315317200000 ;
  	z[5][0] = -.4460530000000E-01 ;
  	u[5][0] = 0.5199212681014;
  	v[5][0] = 0.7198022773915;
  	w[5][0] = 0.2297489167762E-02;

    /* JUPITER */
 	x[6][0] =  4.286062400000;
  	y[6][0] =  -2.621093500000;
  	z[6][0] =  -.8513666500000E-01;
  	u[6][0] =  0.2231109535641;
  	v[6][0] =  0.3949656910948;
  	w[6][0] =  -.6631270424191E-02;

    /* SATURN */
	x[7][0] = 8.834561200000;
  	y[7][0] = 3.097512000000;
  	z[7][0] = -.4051782400000;
  	u[7][0] = -.1244215317121;
  	v[7][0] = 0.3058012696789;
  	w[7][0] = -.3744219597147E-03;

    /* URANUS */
	x[8][0] = 12.29164300000;
  	y[8][0] = -15.57520200000;
  	z[8][0] = -.2172011500000;
  	u[8][0] = 0.1780033999880;
  	v[8][0] = 0.1314573649760;
  	w[8][0] = -.1824027526613E-02;


    /* NEPTUNE */
	x[9][0] = 4.84009700000;
  	y[9][0] = -26.23912700000;
  	z[9][0] = 0.1982557900000;
  	u[9][0] = 0.1578550564045;
  	v[9][0] = 0.9132161165808E-01;
  	w[9][0] = -.5510764371051E-02;

    /* PLUTO */
	x[10][0] = -12.10226300000;
  	y[10][0] = -26.73256000000;
  	z[10][0] = 6.362842900000;
  	u[10][0] = 0.1714028159354;
  	v[10][0] = -.1021868903979;
  	w[10][0] = -.3854629379438E-01;

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

    START_TIMER(main)

    for (int step = 0; step <= num_ts; step++)
    {
        if (debug)
        {
            printf("Step %d\n", step);
        }

        // LOOP 1, you're only assigning once and not touching
	// it again. Parallelizable.
        for (int i = 1; i <= num_bodies; i++)
        {
            // LOOP 2, see above
            for(int j = 1; j <= i - 1; j++)
            {
                X[i][j][0] = x[j][0] - x[i][0];
                Y[i][j][0] = y[j][0] - y[i][0];
                Z[i][j][0] = z[j][0] - z[i][0];
                r[i][j][0] = pow(X[i][j][0], 2) + pow(Y[i][j][0], 2) + pow(Z[i][j][0], 2);
                b[i][j][0] = pow(r[i][j][0], -1.5);
            }

            // LOOP 3, see above
            for (int j = i + 1; j < num_bodies; j++)
            {
                X[i][j][0] = x[j][0] - x[i][0];
                Y[i][j][0] = y[j][0] - y[i][0];
                Z[i][j][0] = z[j][0] - z[i][0];
                r[i][j][0] = pow(X[i][j][0], 2) + pow(Y[i][j][0], 2) + pow(Z[i][j][0], 2);
                b[i][j][0] = pow(r[i][j][0], -1.5);
            }
        }

        // LOOP 4, parallelizable, but a same kind of situation
	// as LOOP 6
        for (int k = 1; k <= mac_degree; k++)
        {
            // LOOP 5 Parallelizable
            for (int i = 1; i <= num_bodies; i++)
            {
                x[i][k] = u[i][k - 1] / k;
                y[i][k] = v[i][k - 1] / k;
                z[i][k] = w[i][k - 1] / k;
            }

            // LOOP 6, parallelizable, but only if we make LOOP 7
            // and LOOP 8 omp single
            for (int i = 1; i <= num_bodies; i++)
            {
                // LOOP 7, not parallelizable, cauchy_prod and cauchy_pow
		// contain loop carry dependencies
                for (int j = 1; j <= num_bodies; j++)
                {
                    X[i][j][k] = x[j][k] - x[i][k];
                    Y[i][j][k] = y[j][k] - y[i][k];
                    Z[i][j][k] = z[j][k] - z[i][k];
                    r[i][j][k] = cauchy_prod(X[i][j], X[i][j], k) + 
                                    cauchy_prod(Y[i][j], Y[i][j], k) + 
                                    cauchy_prod(Z[i][j], Z[i][j], k);
                    b[i][j][k] = cauchy_power(r[i][j], b[i][j], k - 1, -1.5);
                }

                // LOOP 8, see above
                for (int j = i + 1; j <= num_bodies; j++)
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
            // LOOP 9, not parallelizable
            for (int i = 1; i <= num_bodies; i++)
            {
                u[i][k] = 0;
                v[i][k] = 0;
                w[i][k] = 0;
                
                // LOOP 10, see above
                for (int j = 1; j <= i - 1; j++)
                {
                    u[i][k] += mass[j] * cauchy_prod(X[i][j], b[i][j], k - 1);
                    v[i][k] += mass[j] * cauchy_prod(Y[i][j], b[i][j], k - 1);
                    w[i][k] += mass[j] * cauchy_prod(Z[i][j], b[i][j], k - 1);
                }
                
                // LOOP 11, see above
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

        // Determine the values of the Maclaurin polynomial using Horner's algorithm and the
        // stored Maclauren coefficients
	//    
	// Horner_value contains no loop carry dependencies, so this is parallelizable.
        for (int i = 1; i <= num_bodies; i++)
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

        // Output the step number based on the output granularity
        if (step % granularity == 0)
        {   
            if (verbose) 
            {
                printf("Step %d\n", step);
                for (int i = 1; i <= num_bodies; i++)
                {
                    printf("body %d: x: %lf\ty: %lf\tz: %lf\n", i, x[i][0], y[i][0], z[i][0]);
                }
            }
        }
        
        for (int i = 0; i < num_bodies; i++)
        {
            fprintf(body_output[i + 1], "x: %lf\ty: %lf\tz: %lf\n", x[i][0], y[i][0], z[i][0]);
        }
    }
    
    STOP_TIMER(main)
    
    if (timer)
    {
        printf("Main: %0.4fs\n", GET_TIMER(main));
    }
    
    if (debug)
    {
        printf("End of main loop\n\n");
    }
    

    // Cleanup and exit
    if (debug)
    {
        printf("Cleanup and exit\n\n");
    }
    
    for (int i = 0; i < num_bodies; i++)
    {
        fclose(body_output[i + 1]);
    }
    
    return 0;
}