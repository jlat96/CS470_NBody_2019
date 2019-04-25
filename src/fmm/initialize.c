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
    
    int i;
    double mass = 0;
    vector vcm;
    
    double delx, dely, delz;
    
    vcm.x = 0;
    vcm.y = 0;
    vcm.z = 0;
    
    for( i=0; i<N; i++ ){
        
        BD[i].num = i;
        
        BD[i].pos.x = scatter * 0.5 + (drand48()-0.5) * scatter; //not neg so fits planetRegion
        BD[i].pos.y = scatter * 0.5 + (drand48()-0.5) * scatter; //not neg so fits planetRegion
        BD[i].pos.z = scatter * 0.5 + (drand48()-0.5) * scatter; //not neg so fits planetRegion
        
        BD[i].vel.x = v * (drand48()-.5);
        BD[i].vel.y = v * (drand48()-.5);
        BD[i].vel.z = v * (drand48()-.5);
        
        //BD[i].r = 0.3;
        //BD[i].r = 0.03;
        BD[i].r = r_i;
        BD[i].m = 4/3 * M_PI * rho * BD[i].r * BD[i].r * BD[i].r;
        
        mass += BD[i].m;
        vcm   = vcm + BD[i].vel * BD[i].m;
    }
    
    vcm = vcm / mass;
    
    for( i=0; i<N; i++ ){
        BD[i].vel = BD[i].vel - vcm;
    }
}

void InitializeFromFile(planet BD[], char* file_name, double r_i) {
    printf("planet file name: %s\n", file_name);
    FILE *file=fopen(file_name, "r");
    if (file == NULL) {
        printf("Check file name\n");
        exit(EXIT_FAILURE);
    }
    int i = 0;
    double totalMass = 0;
    vector totalVcm; totalVcm.x = 0; totalVcm.y = 0; totalVcm.z = 0;
    double xPos, yPos, zPos, xVel, yVel, zVel, mass;
    while(fscanf(file, "%lf %lf %lf %lf %lf %lf %lf", &xPos, &yPos, &zPos, &xVel, &yVel, &zVel, &mass) != EOF) {
        
        BD[i].num = i;
        
        BD[i].pos.x = xPos;
        BD[i].pos.y = yPos;
        BD[i].pos.z = zPos;
        
        BD[i].vel.x = xVel;
        BD[i].vel.y = yVel;
        BD[i].vel.z = zVel;
        
        BD[i].r = r_i;
        
        BD[i].m = mass;
        
        totalMass += BD[i].m;
        
        totalVcm  += BD[i].vel * BD[i].m;
        
        i++;
    }
    if (i != N) {
        printf("!ERROR: N in config must equal the number of planets in the planet file\n");
        exit(EXIT_FAILURE);
    }
    printf("!Done reading from file\n");
    fclose(file);
    totalVcm = totalVcm / totalMass;
    
    for(int j=0; j<N; j++){
        BD[j].vel = BD[j].vel - totalVcm;
    }
}
