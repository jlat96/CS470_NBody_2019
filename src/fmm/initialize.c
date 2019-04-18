#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "planet.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

extern int N;
extern int N_overflow;
extern double L;
extern double rho;
extern double r_i;

void Initialize( planet BD[], double scatter, double v, double r_i ){
    
    double mass = 0;
    vector vcm;
    
    vcm.x = 0;
    vcm.y = 0;
    vcm.z = 0;
    
    for(int i=0; i<N; i++ ){ // can't be parallelized with omp for because
                             // drand48 relies on previous value, which would
                             // change results in parallel versions of the code
        
        BD[i].num = i;
        
        BD[i].pos.x = scatter * 0.5 + (drand48()-0.5) * scatter;
        BD[i].pos.y = scatter * 0.5 + (drand48()-0.5) * scatter;
        BD[i].pos.z = scatter * 0.5 + (drand48()-0.5) * scatter;
        
        BD[i].vel.x = v * (drand48()-.5);
        BD[i].vel.y = v * (drand48()-.5);
        BD[i].vel.z = v * (drand48()-.5);
        
        BD[i].r = r_i;
        BD[i].m = 4/3 * M_PI * rho * BD[i].r * BD[i].r * BD[i].r;
        
        mass += BD[i].m;
        
        vcm  += BD[i].vel * BD[i].m;
    }
        
    vcm = vcm / mass;

#   pragma omp parallel for default(none) shared(N, BD, vcm)
    for(int i=0; i<N; i++){
        BD[i].vel = BD[i].vel - vcm;
    }
}
