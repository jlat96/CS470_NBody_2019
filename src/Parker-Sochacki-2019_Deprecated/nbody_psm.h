/*
 * nbody_psm.h
 *
 */

#define DEFAULT_TS 0.00390625       // Default time step
#define DEFAULT_NUM 100000          // Defailt number of time steps
#define DEFAULT_GRANULARITY 1000    // Default output granularity
#define MAC_DEGREE 8                // Degree of Maclaurin polyomials

// timing macros (must first declare "struct timeval tv")
#define START_TIMER(NAME) gettimeofday(&tv, NULL); \
    double NAME ## _time = tv.tv_sec+(tv.tv_usec/1000000.0);
#define STOP_TIMER(NAME) gettimeofday(&tv, NULL); \
    NAME ## _time = tv.tv_sec+(tv.tv_usec/1000000.0) - (NAME ## _time);
#define GET_TIMER(NAME) (NAME##_time)

double cauchy_power(double a[], double b[], int n, double pow);
double cauchy_prod(double a[], double b[], int n);
double horner_value(double c[], double t, int n);
void usage(char *argv[]);

