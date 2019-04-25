#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "planet.h"
#include "pln.h"
#include "region.h"
#include "force.h"
#include "exactForce.h"
#include "timing.h"

extern int N;
extern double dt;
extern int collision_check;
extern int collision_number;

/****************************************/
void Reset_Accelerations( planet BD[] ){
    for(int i=0;i<N;i++){
        BD[i].acc.x = 0;
        BD[i].acc.y = 0;
        BD[i].acc.z = 0;
    }
}

/****************************************/
/*        Leapfrog Integrator           */
/****************************************/

/****************************************/

void Position_Half_Step( planet BD[] ){
    
    for(int i=0; i<N; i++ ){
        BD[i].pos = BD[i].pos + BD[i].vel*dt/2;
    }
}

/****************************************/

void Velocity_Full_Step( planet BD[] ){
    
    for(int i=0; i<N; i++ ){
        BD[i].vel = BD[i].vel + BD[i].acc*dt;
    }
}

/****************************************/

void leapfrog( planet BD[], int method )
{
    int i;
    
    Position_Half_Step(BD);
    Reset_Accelerations(BD);
    
    collision_check  = 1;
    collision_number = 0;
    
    if( method == 0 ){
        Exact_Force(BD);
    }
    
    
    else{
        for( i=0; i<N; i++ ){
            //printf("!-----Doing force calculation for planet %d\n",i);
            forceMagic(octree, BD[i], BD);
        }
    }
    collision_check = 0;
    
    Velocity_Full_Step(BD);
    Position_Half_Step(BD);
}


/****************************************/
/*         Omelyan Integrator           */
/****************************************/

double eps = 0.1786178958448;
double lam = -0.2123418310626;
double chi = -0.06626458266982;

//Link to the paper I got this from
//https://arxiv.org/pdf/cond-mat/0110585.pdf

/****************************************/

void Position_Step_1( planet BD[] )
{
    for(int i=0;i<N;i++){
        BD[i].pos = BD[i].pos + BD[i].vel*eps*dt;
    }
}

/****************************************/

void Velocity_Step_1(planet BD[])
{
    for(int i=0;i<N;i++){
        BD[i].vel = BD[i].vel + BD[i].acc*(1-2*lam)*dt/2;
    }
}

/****************************************/

void Position_Step_2(planet BD[])
{
    for(int i=0;i<N;i++){
        BD[i].pos = BD[i].pos + BD[i].vel*chi*dt;
    }
}

/****************************************/

void Velocity_Step_2(planet BD[])
{
    for(int i=0;i<N;i++){
        BD[i].vel = BD[i].vel + BD[i].acc*lam*dt;
    }
}

/****************************************/

void Position_Step_3(planet BD[])
{
    for(int i=0;i<N;i++){
        BD[i].pos = BD[i].pos + BD[i].vel*(1-2*(chi+eps))*dt;
    }
}

/****************************************/

void Velocity_Step_3(planet BD[])
{
    for(int i=0;i<N;i++){
        BD[i].vel = BD[i].vel + BD[i].acc*lam*dt;
    }
}

/****************************************/

void Position_Step_4(planet BD[])
{
    for (int i=0;i<N;i++){
        BD[i].pos = BD[i].pos + BD[i].vel*chi*dt;
    }
}

/****************************************/

void Velocity_Step_Final(planet BD[])
{
    for(int i=0;i<N;i++){
        BD[i].vel = BD[i].vel + BD[i].acc*(1-2*lam)*dt/2;
    }
}

/****************************************/

void Position_Step_Final(planet BD[])
{ 
    for(int i=0;i<N;i++){
        BD[i].pos = BD[i].pos + BD[i].vel*eps*dt;
    }
}

/****************************************/

void omelyan( planet BD[], int method )
{
    //Position_Step_1(BD);
#   pragma omp parallel default(none) shared(BD, N, method, chi, lam, eps, dt, octree, collision_check, collision_number)
    {
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].pos.x = BD[i].pos.x + BD[i].vel.x*eps*dt;
            BD[i].pos.y = BD[i].pos.y + BD[i].vel.y*eps*dt;
            BD[i].pos.z = BD[i].pos.z + BD[i].vel.z*eps*dt;   
        }
        
        //Reset_Accelerations(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].acc.x = 0;
            BD[i].acc.y = 0;
            BD[i].acc.z = 0;
        }
#       pragma omp single
        {
            if( method == 0)
            {
                Exact_Force(BD);
            }
            else
            {
                for(int i=0; i<N; i++)
                {
                    forceMagic(octree, BD[i], BD);
                }
            }
        }
        
        //Velocity_Step_1(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].vel.x = BD[i].vel.x + BD[i].acc.x*(1-2*lam)*dt/2;
            BD[i].vel.y = BD[i].vel.y + BD[i].acc.y*(1-2*lam)*dt/2;
            BD[i].vel.z = BD[i].vel.z + BD[i].acc.z*(1-2*lam)*dt/2;

        }
        
        //Position_Step_2(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].pos.x = BD[i].pos.x + BD[i].vel.x*chi*dt;
            BD[i].pos.y = BD[i].pos.y + BD[i].vel.y*chi*dt;
            BD[i].pos.z = BD[i].pos.z + BD[i].vel.z*chi*dt;

        }
        
        //Reset_Accelerations(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].acc.x = 0;
            BD[i].acc.y = 0;
            BD[i].acc.z = 0;
        }
        
#       pragma omp single
        {
            if( method == 0)
            {
                Exact_Force(BD);
            }
            else
            {
                for(int i=0; i<N; i++)
                {
                    forceMagic(octree, BD[i], BD);
                }
            }
        }
        
        //Velocity_Step_2(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].vel.x = BD[i].vel.x + BD[i].acc.x*lam*dt;
            BD[i].vel.y = BD[i].vel.y + BD[i].acc.y*lam*dt;
            BD[i].vel.z = BD[i].vel.z + BD[i].acc.z*lam*dt;

        }
        
        //Position_Step_3(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].pos.x = BD[i].pos.x + BD[i].vel.x*(1-2*(chi+eps))*dt;
            BD[i].pos.y = BD[i].pos.y + BD[i].vel.y*(1-2*(chi+eps))*dt;
            BD[i].pos.z = BD[i].pos.z + BD[i].vel.z*(1-2*(chi+eps))*dt;

        }
        
        //Reset_Accelerations(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].acc.x = 0;
            BD[i].acc.y = 0;
            BD[i].acc.z = 0;
        }
        
#       pragma omp single
        {
            if( method == 0)
            {
                Exact_Force(BD);
            }
            else
            {
                for(int i=0; i<N; i++)
                {
                    forceMagic(octree, BD[i], BD);
                }
            }
        }
        
        //Velocity_Step_3(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].vel.x = BD[i].vel.x + BD[i].acc.x*lam*dt;
            BD[i].vel.y = BD[i].vel.y + BD[i].acc.y*lam*dt;
            BD[i].vel.z = BD[i].vel.z + BD[i].acc.z*lam*dt;
        }
        
        //Position_Step_4(BD);
#       pragma omp for
        for (int i=0;i<N;i++){
            BD[i].pos.x = BD[i].pos.x + BD[i].vel.x*chi*dt;
            BD[i].pos.y = BD[i].pos.y + BD[i].vel.y*chi*dt;
            BD[i].pos.z = BD[i].pos.z + BD[i].vel.z*chi*dt;

        }
        
        //Reset_Accelerations(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].acc.x = 0;
            BD[i].acc.y = 0;
            BD[i].acc.z = 0;
        }
        
        collision_check  = 1;
        collision_number = 0;
        
#       pragma omp single
        {
            if( method == 0)
            {
                Exact_Force(BD);
            }
            else
            {
                for(int i=0; i<N; i++)
                {
                    forceMagic(octree, BD[i], BD);
                }
            }
        }
        
        collision_check = 0;
        //Velocity_Step_Final(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].vel.x = BD[i].vel.x + BD[i].acc.x*(1-2*lam)*dt/2;
            BD[i].vel.y = BD[i].vel.y + BD[i].acc.y*(1-2*lam)*dt/2;
            BD[i].vel.z = BD[i].vel.z + BD[i].acc.z*(1-2*lam)*dt/2;

        }
        
        //Position_Step_Final(BD);
#       pragma omp for
        for(int i=0;i<N;i++){
            BD[i].pos.x = BD[i].pos.x + BD[i].vel.x*eps*dt;
            BD[i].pos.y = BD[i].pos.y + BD[i].vel.y*eps*dt;
            BD[i].pos.z = BD[i].pos.z + BD[i].vel.z*eps*dt;
        }
    }
}
