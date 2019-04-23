#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "vector.h"
#include "planet.h"
#include "pln.h"
#include "region.h"
#include "collisions.h"
#include "initialize.h"
#include "exactForce.h"
#include "systemProperties.h"
#include "integrators.h"
#include "popRegions.h"
#include "force.h"
#include "resetOctree.h"
#include "timing.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

//UNITS
//Mass   : Earth Mass
//Length : AU (astronomical unit)
//Time   : YEAR
//fix planet sizes, way too big in these units.
//or does this not matter?

int LVL;
int numDel;
int collision_check;
int collision_number;
int collision_pair_1[30000];
int collision_pair_2[30000];
double alpha;
float alphasq;
long int proclist[100];

int N;                   // initial number of planets    
int a;                   // number of inelastic collisions
double dt;               // time step
double r_i;              // initial planet radius
double G  = 4*M_PI*M_PI; // Newton's constant 
double L;                // octree size, dynamically chosen
double rho = 20;         // density of planets

int total_regions;
int currTotalRegions;

vector BDCOM; // center of mass of bodies

region octree;


// calculates center of mass
vector getCoM(planet *BD)
{
    vector com; com.x=0; com.y=0; com.z=0;
    double mTot=0;
#   pragma omp parallel for reduction(+:mTot) default(none) shared(N, com, BD)
    for(int i=0; i<N; i++)
    {
#       pragma omp atomic
        com.x += BD[i].pos.x*BD[i].m;
#       pragma omp atomic
        com.y += BD[i].pos.y*BD[i].m;
#       pragma omp atomic
        com.z += BD[i].pos.z*BD[i].m;
        mTot += BD[i].m;
    }
    com = com / mTot;
    
    return com;
}


void setL(planet *BD, vector com)
{
	//Set L based on particle locations, choosing twice the largest single coordinate relative to center of mass (assuming CoM was properly set to 0,0,0 in initialize.c
	//FIXME what do we do if particle leave region 0? Maybe make a new list separate from the region list to keep track of those?
	//FIXME Should do magnitudes
	double max;
    // does not make sense to parallelize since max could be modified so often
    // that all benefits of parallelization would be cancelled out by constantly
    // locking and unlocking
	for(int i=0; i<N; i++)
	{
        if( fabs(BD[i].pos.x - com.x) > max) max = fabs(BD[i].pos.x - com.x);
        if( fabs(BD[i].pos.y - com.y) > max) max = fabs(BD[i].pos.y - com.y);
        if( fabs(BD[i].pos.z - com.z) > max) max = fabs(BD[i].pos.z - com.z);
    }
	L=max * 4;//  max is the largest single coordinate of a particle. We will center region 0 at 0,0,0 and we want it to be sixteen times as large on every side as the furthest particle
}


void showOctree(region parent)
{
    int i;
    int j=0;
    printf("!P: level=%d, #planets=%d\n", parent.level, parent.numPln);
    for(i=0; i<8; i++){
        if(parent.child[i] != NULL){
            j+=parent.child[i]->numPln;
            printf("!***C: level=%d, #planets=%d\n", parent.child[i]->level, parent.child[i]->numPln);
        }
    }
    printf("!***Diff in planet num = %d\n",j-parent.numPln);
    
    for(i=0;i<8;i++){
        if(parent.child[i] != NULL){
            showOctree(*(parent.child[i]));
        }
    }
}



int main(int nParam, char **paramList)
{  
    long int time_start = clock();
    
    char var[100], val[100];//Placeholders to be used when reading from config.txt
    FILE *config=fopen("config.txt", "r");
    char *file_name = NULL;
    
    double v, scatter, r_i, t_end;
    int debug;
    int fSkip;
    
    while( fscanf(config, "%s %s", var, val) != EOF)
    {
        if( strcmp(var, "LVL") == 0)       LVL       = atoi(val);
        if( strcmp(var, "N") == 0)         N         = atoi(val);
        if( strcmp(var, "dt") == 0)        dt        = atof(val);
        if( strcmp(var, "frameskip") == 0) fSkip     = atoi(val);
        if( strcmp(var, "v_i") == 0)       v         = atof(val);
        if( strcmp(var, "scatter") == 0)   scatter   = atof(val);
        if( strcmp(var, "alpha") == 0)     alpha     = atof(val);
        if( strcmp(var, "r_i") == 0)       r_i       = atof(val);
        if( strcmp(var, "t_end") == 0)     t_end     = atof(val);
        if( strcmp(var, "debug") == 0)     debug     = atof(val); 
        if( strcmp(var, "file_name") == 0) file_name =      val ;
    }
    
    int fil=1;
    while(fil<nParam-1)
    {
        if (debug) {
            printf("!Reading command line params  ..  fil:%d  pL:%s  pL+1:%s \n", fil, paramList[fil], paramList[fil+1]);
        }
        if( strcmp(paramList[fil], "LVL") == 0)       LVL       = atoi(paramList[fil+1]);
        if( strcmp(paramList[fil], "N") == 0)         N         = atoi(paramList[fil+1]);
        if( strcmp(paramList[fil], "dt") == 0)        dt        = atof(paramList[fil+1]);
        if( strcmp(paramList[fil], "frameskip") == 0) fSkip     = atoi(paramList[fil+1]);
        if( strcmp(paramList[fil], "v_i") == 0)       v         = atof(paramList[fil+1]);
        if( strcmp(paramList[fil], "scatter") == 0)   scatter   = atof(paramList[fil+1]);
        if( strcmp(paramList[fil], "alpha") == 0)     alpha     = atof(paramList[fil+1]);
        if( strcmp(paramList[fil], "r_i") == 0)       r_i       = atof(paramList[fil+1]);
        if( strcmp(paramList[fil], "t_end") == 0)     t_end     = atof(paramList[fil+1]);
        if( strcmp(paramList[fil], "debug") == 0)     debug     = atof(paramList[fil+1]);  
        if( strcmp(paramList[fil], "file_name") == 0) file_name =      paramList[fil+1] ;
        fil+=2;
    }
    
    planet BD[N];
    
    int i;
    
    int frame    = 0;
    int method   = 1;
    
    double t;
    double M_init, E_init, P_init;
    double M_curr, E_curr, P_curr;
    
    if (debug) {
        if(method == 1){ printf("!method = 1. Using multipole method\n");}
        else if(method == 0){ printf("!method = 0. Using exact force\n");} 
        else {printf("!Invalid entry for force method. Please choose either 0 or 1\n"); exit(0);}
    }
    alphasq = alpha*alpha;
    
    for(i=0;i<100;i++) proclist[i]=0;
    
    if (file_name != NULL) {
        InitializeFromFile(BD, file_name ,r_i);
    } else {
        //TODO add timer code and compare to serial version
        Initialize(BD, scatter, v, r_i);
    }
    
    if (debug) {
        printf("!Done intializing planets\n");
    }
    
    if (debug) {
        printf("!Calculating center of mass\n");
    }
    //TODO add timer code and compare to serial version
    BDCOM = getCoM(BD);
    
    setL(BD, BDCOM);
    if (debug) {
        printf("!fmm.c: L = %1.3e\n", L);
    }
    
    octree.numPln = 0;
    octree.level  = 0;
    octree.size   = L;
    
    octree.location.x = BDCOM.x-L/2.0;
    octree.location.y = BDCOM.y-L/2.0;
    octree.location.z = BDCOM.z-L/2.0;
    
        
    //TODO add timer code and compare to serial version
    M_init = Mass(BD);
    //TODO add timer code and compare to serial version
    E_init = Energy(BD);
    //TODO add timer code and compare to serial version
    P_init = Momentum(BD);
    
    if (debug) {
        printf("Initial mass: %f\n", M_init);
        printf("Initial energy: %f\n", E_init);
        printf("Initial momentum: %f\n", P_init);
    }
    
    M_curr = M_init;
    E_curr = E_init;
    P_curr = P_init;
    
    if (debug) {
        printf("!fmm.c: Done everything before the time loop\n");
    }
    for( t=0; t <= t_end; t+=dt )
    {  
        
        if(t!=0) {
            //TODO add timer code and compare to serial version
            resetOctree(octree);
        }
        pop_level_0(octree, BD);
        
        total_regions=1;
        
        loopOverRegions(octree, BD);
        
        //TODO add timer code and compare to serial version
        recurse_divide_by_mass(octree);
        
        //TODO add timer code and compare to serial version
        omelyan(BD, method);
        
        //TODO add timer code and compare to serial version
        a += collisionCheck(BD);
        
        if( frame % fSkip == 0 )
        {
            if( frame % 1 == 0 )
            {
                //TODO add timer code and compare to serial version
                M_curr = Mass(BD);
                //TODO add timer code and compare to serial version
                E_curr = Energy(BD);
                //TODO add timer code and compare to serial version
                P_curr = Momentum(BD);
            }
                        
            if (debug) {
                printf("E(0) = %.2f : E(t) = %.2f : P = %.2f : t = %.2f : Planets = %d : Mass = %.2f\n", E_init, E_curr, P_curr, t, N, M_curr);
                for( i=0; i<N; i++ ) {
                    printf("c3 %e %e %e %e\n", BD[i].pos.x, BD[i].pos.y, BD[i].pos.z, BD[i].r);
                }
                printf("F\n");   
            }
        }
        frame++;
    }
    printf("E(0) = %.2f : E(t) = %.2f : P = %.2f : t = %.2f : Planets = %d : Mass = %.2f\n", E_init, E_curr, P_curr, t_end, N, M_curr);
    printf("time: %f seconds\n", (clock() - time_start) / (double) CLOCKS_PER_SEC);
}
