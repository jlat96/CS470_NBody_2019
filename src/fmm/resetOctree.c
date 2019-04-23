#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pln.h"
#include "vector.h"
#include "planet.h"
#include "region.h"
#include "popRegions.h"

#ifdef _OPENMP
    #include <omp.h>
#endif
extern int LVL;

void resetOctree(region& parent_region){
    //printf("Parent is on level %d\n", parent_region.level);
    if( parent_region.planets != NULL ) clearList(&parent_region.planets);
    
    //if(parent_region.child[0] != NULL && parent_region.level != LVL-1){
#   pragma omp parallel for default(none) shared(parent_region)
    for(int i=0; i<8; i++){
        if(parent_region.child[i]!=NULL)
        { // it has kids to erase
            resetOctree(*parent_region.child[i]);
            free(parent_region.child[i]);
            parent_region.child[i]=NULL;
        }
    }

}

