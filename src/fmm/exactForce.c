#include "vector.h"
#include "planet.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

extern int N;
extern double G;

void Exact_Force( planet BD[] ){
    
#   pragma omp parallel for default(none) shared(N, BD, G)
    for(int i=0; i<N; i++ ){
        for(int j=i+1; j<N; j++){
            double rad = radius( BD[i].pos, BD[j].pos );
            vector F   = G * BD[j].m * BD[i].m * (BD[j].pos-BD[i].pos) / (rad*rad*rad);
            
#           pragma omp atomic
            BD[i].acc.x = BD[i].acc.x + F.x/BD[i].m;
#           pragma omp atomic
            BD[i].acc.y = BD[i].acc.y + F.y/BD[i].m;
#           pragma omp atomic
            BD[i].acc.z = BD[i].acc.z + F.z/BD[i].m;
#           pragma omp atomic
            BD[j].acc.x = BD[j].acc.x + F.x/BD[j].m;
#           pragma omp atomic
            BD[j].acc.y = BD[j].acc.y + F.y/BD[j].m;
#           pragma omp atomic
            BD[j].acc.z = BD[j].acc.z + F.z/BD[j].m;
        }
    }
}
