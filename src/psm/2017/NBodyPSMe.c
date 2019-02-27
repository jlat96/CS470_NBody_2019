/*
 This program is a Picard-Maclaurin solver for Newton's N Body Problem.
 It is currently set up for the 10 planet solar system.
 The Maclaurin polynomial are determined from Picard iteration using Cauchy
 products. The algorithm is at http://educ.jmu.edu/~sochacjs/PSM.html -
 An Expository Document on Using the Modified Picard Method to Solve Initial
 Value Ordinary Differential and Partial Differential Equations
 and Planetary Motion Application Using Parker Sochacki
 
 The input for this code is the number of bodies (NB), the mass of each body,
 the initial position (x[i,0],y[i,0],z[i,0]) and the initial velocity
 (u[i,0],v[i,0],w[i,0]) of each body, the time step,  the time to
 run the simulation and the degree of the Maclaurin polynomial for the
 numerical solution. There is a special routine for raising a power series to an
 exponent (cauchypower in  CauchyProdPower.h). This reduces the number of
 calculations needed.
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "NBPSMparameters.h"      // Power Series Method parameters

#include "Horner.h"       // Evaluate the Maclaurin Polynomials using Horner's algorithm

#include "CauchyProd.h"   // Cauchy Product routine
#include "CauchyProdPower.h"   // Cauchy Product Power routine

int main()
{
// NB is the numbe of bodies (planets)
	int i,j,k,m,ns,nc;
	double x[NB+1][Mdeg+1], y[NB+1][Mdeg+1], z[NB+1][Mdeg+1];  // position of bodies
	double u[NB+1][Mdeg+1], v[NB+1][Mdeg+1], w[NB+1][Mdeg+1];  // velocity of bodies
    double mass[NB+1];  // mass of the body
    double X[NB+1][NB+1][Mdeg+1], Y[NB+1][NB+1][Mdeg+1], Z[NB+1][NB+1][Mdeg+1];
    double r[NB+1][NB+1][Mdeg+1],b[NB+1][NB+1][Mdeg+1];  // used to help in calculating the body data
    double xPSM,yPSM,zPSM,uPSM,vPSM,wPSM;  // used to update the body location and velocity
    
	FILE *fpx[NB+1],*fpy[NB+1],*fpz[NB+1];
    char filename[10];

    printf(" h = %f \n",h);
    printf(" num = %d \n",num);
    
	/*  The masses. */

/*  mass of body Sun        =  1.
  mass of body Mercury    =  1.66013679527193035E-7
  mass of body Venus      =  2.44783959796682464E-6
  mass of body Earth      =  3.04043273871083524E-6
  mass of body Mars       =  3.22714936215392876E-7
  mass of body Jupiter    =  9.54790662147324233E-4
  mass of body Saturn     =  2.85877644368210402E-4
  mass of body Uranus     =  4.35540069686411149E-5
  mass of body Neptune    =  5.17759138448793649E-5
  mass of body Pluto      =  7.6923076923076926E-9 */

  mass[1]       =  1.;
  mass[2]    =  1.66013679527193035E-7;
  mass[3]      =  2.44783959796682464E-6;
  mass[4]     =  3.04043273871083524E-6;
  mass[5]      =  3.22714936215392876E-7;
  mass[6]    =  9.54790662147324233E-4;
  mass[7]     =  2.85877644368210402E-4;
  mass[8]     =  4.35540069686411149E-5;
  mass[9]    =  5.17759138448793649E-5;
  mass[10]      =  7.6923076923076926E-9;

	/*    Set up the initial positions and velocities  */

/* SUN */

    x[1][0] = 0;
  	y[1][0] = 0;
  	z[1][0] = 0;
  	u[1][0] = 0;
  	v[1][0] = 0;
  	w[1][0] = 0;

/* MERCURY */

  	x[2][0] =  0.1633519000000E-02;
  	y[2][0] =  0.3077199200000;
  	z[2][0] =  0.2498653800000E-01;
  	u[2][0] = -1.963534107001;
  	v[2][0] =  0.6841909835729E-01;
  	w[2][0] = 0.1858224458145;

/* VENUS */

	x[3][0] =  0.2640622900000;
  	y[3][0] = 0.6709840300000;
  	z[3][0] =  -.6078524800000E-02;
  	u[3][0] =  -1.097989091627;
  	v[3][0] =  0.4250727728824;
  	w[3][0] = 0.6918653377540E-01;
	
/* EARTH */

	x[4][0] =  0.6643154100000E-01;
  	y[4][0] =    0.9817401900000;
  	z[4][0] =  0.6625301000000E-05;
  	u[4][0] =  -1.014141067952;
  	v[4][0] =  0.6377084582524E-01;
  	w[4][0] =  -.1047343992878E-05;

/* MARS */ 

	x[5][0] =   1.105892500000;
  	y[5][0] =  -.8315317200000 ;
  	z[5][0] = -.4460530000000E-01 ;
  	u[5][0] = 0.5199212681014;
  	v[5][0] =  0.7198022773915;
  	w[5][0] =  0.2297489167762E-02;

/* JUPITER */

 	x[6][0] =  4.286062400000;
  	y[6][0] =  -2.621093500000;
  	z[6][0] =  -.8513666500000E-01;
  	u[6][0] =  0.2231109535641;
  	v[6][0] =  0.3949656910948;
  	w[6][0] =  -.6631270424191E-02;

/* SATURN */

	x[7][0] =  8.834561200000;
  	y[7][0] =  3.097512000000;
  	z[7][0] = -.4051782400000;
  	u[7][0] = -.1244215317121;
  	v[7][0] =  0.3058012696789;
  	w[7][0] = -.3744219597147E-03;

/* URANUS */

	x[8][0] =  12.29164300000;
  	y[8][0] = -15.57520200000;
  	z[8][0] = -.2172011500000;
  	u[8][0] =  0.1780033999880;
  	v[8][0] =  0.1314573649760;
  	w[8][0] = -.1824027526613E-02;


/* NEPTUNE */

	x[9][0] =  4.84009700000;
  	y[9][0] = -26.23912700000;
  	z[9][0] = 0.1982557900000;
  	u[9][0] =  0.1578550564045;
  	v[9][0] =  0.9132161165808E-01;
  	w[9][0] =  -.5510764371051E-02;

/* PLUTO */

	x[10][0] =  -12.10226300000;
  	y[10][0] = -26.73256000000;
  	z[10][0] = 6.362842900000;
  	u[10][0] =  0.1714028159354;
  	v[10][0] =  -.1021868903979;
  	w[10][0] = -.3854629379438E-01;

   	/* Put initial conditions in output files */

  	for(i=1; i<=NB; i++)
  	{
    		sprintf(filename,"x%d",i);
        	if ((fpx[i] = fopen(filename, "w")) == NULL)
        	{
            		printf("%s not opened\n", filename);
            		exit(1);
        	}
        	fprintf(fpx[i],"%lf\n", x[i][0]);

		sprintf(filename,"y%d",i);
        	if ((fpy[i] = fopen(filename, "w")) == NULL)
         	{
            		printf("%s not opened\n", filename);
            		exit(1);
         	}
         	fprintf(fpy[i],"%lf\n", y[i][0]);

		sprintf(filename,"z%d",i);
         	if ((fpz[i] = fopen(filename, "w")) == NULL)
         	{
            		printf("%s not opened\n", filename);
            		exit(1);
         	}
         	fprintf(fpz[i],"%lf\n", z[i][0]);
  	}

   	/*  Do the  simulation for num time steps.  */

  nc = 0;

  for ( ns = 2; ns <= num; ns++ )
  {
      
      /*
       
       Set up X,Y,Z, r and b from the 'updated' initial conditions
       These are needed for the gravitational force laws (Newton)
       
       */
      
      for ( i = 1; i <=NB; i++ )
      {
          for ( j = 1; j <= i-1; j++ )
          {
              
              X[i][j][0] = x[j][0]-x[i][0];
              Y[i][j][0] = y[j][0]-y[i][0];
              Z[i][j][0] = z[j][0]-z[i][0];
              r[i][j][0] = pow(X[i][j][0],2)+pow(Y[i][j][0],2)+pow(Z[i][j][0],2);
              b[i][j][0] = pow(r[i][j][0],-1.5);
              
          }
          for ( j = i+1; j <= NB; j++ )
          {
              
              X[i][j][0] = x[j][0]-x[i][0];
              Y[i][j][0] = y[j][0]-y[i][0];
              Z[i][j][0] = z[j][0]-z[i][0];
              r[i][j][0] = pow(X[i][j][0],2)+pow(Y[i][j][0],2)+pow(Z[i][j][0],2);
              b[i][j][0] = pow(r[i][j][0],-1.5);
              
          }
      }
      
    for ( k = 1; k <= Mdeg; k++ )
    {
        for ( i = 1; i <= NB; i++ )
        {
           x[i][k] = u[i][k-1]/k;
           y[i][k] = v[i][k-1]/k;
           z[i][k] = w[i][k-1]/k;
        }
        for ( i = 1; i <= NB; i++ )
        {
           for ( j = 1; j <= i-1; j++ )
           {
                
                X[i][j][k] = x[j][k]-x[i][k];
                Y[i][j][k] = y[j][k]-y[i][k];
                Z[i][j][k] = z[j][k]-z[i][k];
                r[i][j][k] = cauchyprod(X[i][j],X[i][j],k)+cauchyprod(Y[i][j],Y[i][j],k)+cauchyprod(Z[i][j],Z[i][j],k);
                b[i][j][k] = cauchypower(r[i][j],b[i][j],k-1,-1.5);
                
            }
            for ( j = i+1; j <= NB; j++ )
            {
                
                X[i][j][k] = x[j][k]-x[i][k];
                Y[i][j][k] = y[j][k]-y[i][k];
                Z[i][j][k] = z[j][k]-z[i][k];
                r[i][j][k] = cauchyprod(X[i][j],X[i][j],k)+cauchyprod(Y[i][j],Y[i][j],k)+cauchyprod(Z[i][j],Z[i][j],k);
                b[i][j][k] = cauchypower(r[i][j],b[i][j],k-1,-1.5);
                
            }
        }
        for ( i = 1; i <= NB; i++ )
        {
            u[i][k] = 0;
            v[i][k] = 0;
            w[i][k] = 0;
            
            for ( j = 1; j <= i-1; j++ )
            {
                
                u[i][k] = u[i][k] + mass[j]*cauchyprod(X[i][j],b[i][j],k-1);
                v[i][k] = v[i][k] + mass[j]*cauchyprod(Y[i][j],b[i][j],k-1);
                w[i][k] = w[i][k] + mass[j]*cauchyprod(Z[i][j],b[i][j],k-1);
                
            }
            for ( j = i+1; j <= NB; j++ )
            {
                
                u[i][k] = u[i][k] + mass[j]*cauchyprod(X[i][j],b[i][j],k-1);
                v[i][k] = v[i][k] + mass[j]*cauchyprod(Y[i][j],b[i][j],k-1);
                w[i][k] = w[i][k] + mass[j]*cauchyprod(Z[i][j],b[i][j],k-1);
                
            }
            u[i][k] = u[i][k]/k;
            v[i][k] = v[i][k]/k;
            w[i][k] = w[i][k]/k;
      
        }
        
      }
      
      // Get the values of the Maclaurin polynomials using Horner's algorithm
      
      for ( i = 1; i <= NB; i++ )
      {
          xPSM=horner(x[i],h,Mdeg);  // x coordinate of body
          x[i][0]=xPSM;
          yPSM=horner(y[i],h,Mdeg);  // y coordinate of body
          y[i][0]=yPSM;
          zPSM=horner(z[i],h,Mdeg);  // z coordinate of body
          z[i][0]=zPSM;
          uPSM=horner(u[i],h,Mdeg);  // x velociy of body
          u[i][0]=uPSM;
          vPSM=horner(v[i],h,Mdeg);  // y velociy of body
          v[i][0]=vPSM;
          wPSM=horner(w[i],h,Mdeg);  // z velociy of body
          w[i][0]=wPSM;
      }
        
	nc = nc+1;
      
	if (nc==np)
	{
	  printf(" step = %d \n",ns);
	  nc = 0;
    }
    for(i=1; i<=NB; i++)	/* output the (x,y,z) position of each body */
    {

        fprintf(fpx[i],"%lf\n",x[i][0]);
        fprintf(fpy[i],"%lf\n",y[i][0]);
        fprintf(fpz[i],"%lf\n",z[i][0]);

    }

  }


  	for(i=1; i<=NB; i++)	/* close the position files */
	{
      	fclose(fpx[i]);
        fclose(fpy[i]);
        fclose(fpz[i]);
    }

    return(0);

}

