/*
 * nbody_psm.h
 *
 */

#define DEFAULT_TS pow(2, -8)       // Default time step
#define DEFAULT_NUM 100000          // Defailt number of time steps
#define DEFAULT_GRANULARITY 1000    // Default output granularity
#define MAC_DEGREE 8                // Degree of Maclaurin polyomials
#define NUM_BODIES 10               // Number of Bodies

double cauchy_power(double a[], double b[], int n, double pow);
double cauchy_prod(double a[], double b[], int n);
double cauchy_square(double a[], int n);
double horner_value(double c[], double t, int n);
void usage(char *argv[]);

