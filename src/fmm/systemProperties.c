#include "vector.h"
#include "planet.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

extern int N;
extern double G;

/**********************************************/

double Energy(planet BD[]){
    
    double KE = 0, PE = 0;
    
#   pragma omp parallel default(none) shared(KE, PE, N, BD)
    {
#       pragma omp for reduction(+:KE)
        for(int i=0; i<N; i++ ){
            KE += 0.5 * BD[i].m * ( BD[i].vel * BD[i].vel );
        }

#       pragma omp for reduction(+:PE)
        for(int i=0; i<N; i++){
            for(int j=i+1; j<N; j++ ){
                PE += BD[i].m * BD[j].m / radius( BD[i].pos, BD[j].pos );
            }
        }
    }
    return KE - G*PE;
}

/*********************************************/

double Momentum( planet BD[] ){
    
    vector p;
    p.x = 0;
    p.y = 0;
    p.z = 0;
    
#   pragma omp parallel for default(none) shared(p, BD, N)
    for(int i=0; i<N; i++ ){
#       pragma omp atomic
        p.x = p.x + BD[i].m * BD[i].vel.x;
#       pragma omp atomic
        p.y = p.y + BD[i].m * BD[i].vel.y;
#       pragma omp atomic
        p.z = p.z + BD[i].m * BD[i].vel.z;
    }
    
    return Magnitude( p );
}

/********************************************/

double Mass( planet BD[] ){
    
    double mass = 0.0;
#   pragma omp parallel for default(none) reduction(+:mass) shared(N, BD)
    for(int i=0; i<N; i++ ){
        mass += BD[i].m;
    }
    return mass;
}
