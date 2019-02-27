function [x, y, z] = NBodyPSM();

%   This a modified Picard iterator that generates Maclaurin polynomial
%   approximations to the solution of Newton's N body problem under 
%   gravitation. It is set up for three bodies, but one can change it
%   to model the motion of planets by using the data below

%   N below is the number of bodies
%   x(j,k)  j = 1,2, ... , N is the x coordinates of body j
%   y(j,k)  j = 1,2, ... , N is the y coordinates of body j
%   z(j,k)  j = 1,2, ... , N is the z coordinates of body j

% below you enter the information to generate these coordinate values

% If you want to do  the solar system this is the data
  

%{
  mass of body Sun        =  1.
  mass of body Mercury    =  1.66013679527193035E-7
  mass of body Venus      =  2.44783959796682464E-6
  mass of body Earth      =  3.04043273871083524E-6
  mass of body Mars       =  3.22714936215392876E-7
  mass of body Jupiter    =  9.54790662147324233E-4
  mass of body Saturn     =  2.85877644368210402E-4
  mass of body Uranus     =  4.35540069686411149E-5
  mass of body Neptune    =  5.17759138448793649E-5
  mass of body Pluto      =  7.6923076923076926E-9 */


  mass(1)       =  1.;
  mass(2)    =  1.66013679527193035E-7;
  mass(3)      =  2.44783959796682464E-6;
  mass(4)     =  3.04043273871083524E-6;
  mass(5)      =  3.22714936215392876E-7;
  mass(6)    =  9.54790662147324233E-4;
  mass(7)     =  2.85877644368210402E-4;
  mass(8)     =  4.35540069686411149E-5;
  mass(9)    =  5.17759138448793649E-5;
  mass(10)      =  7.6923076923076926E-9;

	%  /*    Set up the initial conditions  */

% /* SUN */

        alpha(1,1) = 0;
  	beta(1,1) = 0;
  	Gamma(1,1) = 0;
  	delta(1,1) = 0;
  	rho(1,1) = 0;
  	lambda(1,1) = 0;

% /* MERCURY */

  	alpha(2,1) =  0.1633519000000E-02;
  	beta(2,1) =  0.3077199200000;
  	Gamma(2,1) =  0.2498653800000E-01;
  	delta(2,1) = -1.963534107001;
  	rho(2,1) =  0.6841909835729E-01;
  	lambda(2,1) = 0.1858224458145;

% /* VENUS */

	alpha(3,1) =  0.2640622900000;
  	beta(3,1) = 0.6709840300000;
  	Gamma(3,1) =  -.6078524800000E-02;
  	delta(3,1) =  -1.097989091627;
  	rho(3,1) =  0.4250727728824;
  	lambda(3,1) = 0.6918653377540E-01;
	
% /* EARTH */

	alpha(4,1) =  0.6643154100000E-01; 
  	beta(4,1) =    0.9817401900000;
  	Gamma(4,1) =  0.6625301000000E-05;
  	delta(4,1) =  -1.014141067952;
  	rho(4,1) =  0.6377084582524E-01;
  	lambda(4,1) =  -.1047343992878E-05;

% /* MARS */ 

	alpha(5,1) =     1.105892500000;
  	beta(5,1) =  -.8315317200000 ;
  	Gamma(5,1) = -.4460530000000E-01 ;
  	delta(5,1) = 0.5199212681014;
  	rho(5,1) =  0.7198022773915;
  	lambda(5,1) =  0.2297489167762E-02;

% /* JUPITER */

 	alpha(6,1) =  4.286062400000; 
  	beta(6,1) =  -2.621093500000;
  	Gamma(6,1) =  -.8513666500000E-01;
  	delta(6,1) =  0.2231109535641;
  	rho(6,1) =  0.3949656910948;
  	lambda(6,1) =  -.6631270424191E-02;

% /* SATURN */

	alpha(7,1) =  8.834561200000;
  	beta(7,1) =  3.097512000000;
  	Gamma(7,1) = -.4051782400000;
  	delta(7,1) = -.1244215317121;
  	rho(7,1) =  0.3058012696789;
  	lambda(7,1) = -.3744219597147E-03;

% /* URANUS */

	alpha(8,1) =  12.29164300000;
  	beta(8,1) = -15.57520200000;
  	Gamma(8,1) = -.2172011500000;
  	delta(8,1) =  0.1780033999880;
  	rho(8,1) =  0.1314573649760;
  	lambda(8,1) = -.1824027526613E-02;


% /* NEPTUNE */

	alpha(9,1) =  4.84009700000;
  	beta(9,1) = -26.23912700000;
  	Gamma(9,1) = 0.1982557900000;
  	delta(9,1) =  0.1578550564045;
  	rho(9,1) =  0.9132161165808E-01;
  	lambda(9,1) =  -.5510764371051E-02;

% /* PLUTO */

	alpha(10,1) =  -12.10226300000;
  	beta(10,1) = -26.73256000000;
  	Gamma(10,1) = 6.362842900000;
  	delta(10,1) =  0.1714028159354;
  	rho(10,1) =  -.1021868903979;
  	lambda(10,1) = -.3854629379438E-01;
%}
  
  N = 3;
  

%    numeq is the number of equations

  numeq = N*6+N*(N-1)/2;
  
%   K1 is the degree of the Maclaurin polynomials that will be used to approximate each
%   component of the solution.

  K1 = 4;
  
%  Give the time step

  h = 1./2^3;
  
%   Give the time

  tyme = 250.;
  
%  The number of time steps can be calculated or chosen

   num_time_steps = ceil(tyme/h)+1;
%  num_time_steps = 100;
  
% np is used for printing out the current time step of the iteration

  np = 100;

%  Give the masses of each body.
  
  mass(1) = 999;
  mass(2) = 1;
  mass(3) = 1;
   
%  Set up the initial conditions

%  alpha(i,.) is the x coordinate of body i   . is for the degree of the Maclaurin polynomial
%  beta(i,.) is the y coordinate of body i   . is for the degree of the Maclaurin polynomial
%  Gamma(i,.) is the z coordinate of body i   . is for the degree of the Maclaurin polynomial
%  delta(i,.) is the x component of the velocity of body i   
%  rho(i,.) is the y component of the velocity of body i   
%  lambda(i,.) is the z component of the velocity of body i   
%  These must be given.
 
  alpha(1,1) = 0;
  beta(1,1) = 0;
  Gamma(1,1) = 0;
  delta(1,1) = 0;
  rho(1,1) = 0;
  lambda(1,1) = 0;
  alpha(2,1) = 100.1;
  beta(2,1) = 0;
  Gamma(2,1) = 0;
  delta(2,1) = 0;
  rho(2,1) = -1;
  lambda(2,1) = 1;
  alpha(3,1) = 0;
  beta(3,1) = 150;
  Gamma(3,1) = 100;
  delta(3,1) = -1;
  rho(3,1) = 0;
  lambda(3,1) = 1;
  
  % the input data is entered. The code generates the coordinates
  % through the algorithm below. Running this code will give you
  % the vectors x,y,z for each body  j = 1,2, ... , N.
  
  % sigma is the symmetric component of the N Body Problem. See www.math.jmu.edu/~jim/picard.htm
  
  for i = 1 : N 
  for j = i+1 : N 
    sigma(i,j,1)= ((alpha(i,1)-alpha(j,1))^2+(beta(i,1)-beta(j,1))^2+(Gamma(i,1)-Gamma(j,1))^2)^(-1/2.);
  end
  end
  
  for i = 1 : N 
  for j = i+1 : N 
    sigma(j,i,1) = sigma(i,j,1);
  end
  end
  
%  Save the positions of each body for plotting below.  
  
  for j = 1 : N
  
    x(j,1) = alpha(j,1);
    y(j,1) = beta(j,1);
    z(j,1) = Gamma(j,1);
  
  end
  
   
%  Do the  simulation.  Here is where the coefficients of the Maclaurin polynomials for the
%  approximations are determined. See www.math.jmu.edu/~jim/picard.htm for the algorithm.

  for ns = 2 : num_time_steps 
    if (mod(ns,np)==0)
        ns
    end
    for k = 1 : K1 
      for i = 1 : N     
        alpha(i,k+1) = delta(i,k)/k;
        beta(i,k+1) = rho(i,k)/k;
        Gamma(i,k+1) = lambda(i,k)/k;
        for j = i+1 : N 
          F(i,j,k) = 0;
          for m = 0 : k-1 
            F(i,j,k) = F(i,j,k) + sigma(i,j,m+1)*sigma(i,j,k-m);
          end
        end
        for j = 1 : i-1 
          F(i,j,k) = F(j,i,k);
        end
        for j = i+1 : N 
          G(i,j,k) = 0;
          for m = 0 : k-1 
            G(i,j,k) = G(i,j,k) + sigma(i,j,m+1)*F(i,j,k-m);
          end
        end
        for j = 1 : i-1 
          G(i,j,k) = G(j,i,k);
        end
        delta(i,k+1) = 0;
        for j = 1 : i-1 
          Q = 0;
          for m = 0 : k-1 
            Q = Q + (alpha(i,m+1)-alpha(j,m+1))*G(i,j,k-m);
          end
          delta(i,k+1) = delta(i,k+1) - mass(j)*Q;
        end
        for j = i+1 : N 
          Q = 0;
          for m = 0 : k-1 
            Q = Q + (alpha(i,m+1)-alpha(j,m+1))*G(i,j,k-m);
          end      
          delta(i,k+1) = delta(i,k+1) - mass(j)*Q;
        end
        delta(i,k+1) = delta(i,k+1)/k;
        rho(i,k+1) = 0;
        for j = 1 : i-1 
          Q = 0;
          for m = 0 : k-1 
            Q = Q + (beta(i,m+1)-beta(j,m+1))*G(i,j,k-m);
          end
          rho(i,k+1) = rho(i,k+1) - mass(j)*Q;
        end
        for j = i+1 : N 
          Q = 0;
          for m = 0 : k-1 
            Q = Q + (beta(i,m+1)-beta(j,m+1))*G(i,j,k-m);
          end      
          rho(i,k+1) = rho(i,k+1) - mass(j)*Q;
        end
        rho(i,k+1) = rho(i,k+1)/k;
        lambda(i,k+1) = 0;
        for j = 1 : i-1 
          Q = 0;
          for m = 0 : k-1 
            Q = Q + (Gamma(i,m+1)-Gamma(j,m+1))*G(i,j,k-m);
          end
          lambda(i,k+1) = lambda(i,k+1) - mass(j)*Q;
        end
        for j = i+1 : N 
          Q = 0;
          for m = 0 : k-1 
            Q = Q + (Gamma(i,m+1)-Gamma(j,m+1))*G(i,j,k-m);
          end      
          lambda(i,k+1) = lambda(i,k+1) - mass(j)*Q;
        end
        lambda(i,k+1) = lambda(i,k+1)/k;
        for j = i+1 : N 
          H(i,j,k) = 0;
          for m = 0 : k-1 
            H(i,j,k) = H(i,j,k) + (alpha(i,m+1)-alpha(j,m+1))*(delta(i,k-m)-delta(j,k-m)) + (beta(i,m+1)-beta(j,m+1))*(rho(i,k-m)-rho(j,k-m)) + (Gamma(i,m+1)-Gamma(j,m+1))*(lambda(i,k-m)-lambda(j,k-m));
          end
          T = 0;
          for m = 0 : k-1 
            T = T + H(i,j,m+1)*G(i,j,k-m);
          end
          sigma(i,j,k+1) = -T/k;
        end      
        for j = 1 : i-1 
          sigma(i,j,k+1) = sigma(j,i,k+1);
        end
      end
    end
    
% Use Horner's algorithm to determine the values of the Maclaurin polynomials using the
% Maclaurin coefficients calculated above.

    for i = 1 : N 
      temp = alpha(i,K1)+alpha(i,K1+1)*h;
      for k = 1 : K1-1 
        temp = temp*h+alpha(i,K1-k);
      end
      alpha(i,1) = temp;
      temp = beta(i,K1)+beta(i,K1+1)*h;
      for k = 1 : K1-1 
        temp = temp*h+beta(i,K1-k);
      end
      beta(i,1) = temp;
      temp = Gamma(i,K1)+Gamma(i,K1+1)*h;
      for k = 1 : K1-1 
        temp = temp*h+Gamma(i,K1-k);
      end
      Gamma(i,1) = temp;
      temp = delta(i,K1)+delta(i,K1+1)*h;
      for k = 1 : K1-1 
        temp = temp*h+delta(i,K1-k);
      end
      delta(i,1) = temp;
      temp = rho(i,K1)+rho(i,K1+1)*h;
      for k = 1 : K1-1 
        temp = temp*h+rho(i,K1-k);
      end
      rho(i,1) = temp;
      temp = lambda(i,K1)+lambda(i,K1+1)*h;
      for k = 1 : K1-1 
        temp = temp*h+lambda(i,K1-k);
      end
      lambda(i,1) = temp;
      for j = i+1 : N 
        temp = sigma(i,j,K1)+sigma(i,j,K1+1)*h;
        for k = 1 : K1-1 
          temp = temp*h+sigma(i,j,K1-k);
        end
        sigma(i,j,1) = temp;
      end
    end
    for i = 1 : N 
    for j = 1 : i-1 
      sigma(i,j,1) = sigma(j,i,1);
    end
    end

%  Store the positions of each body.

     for j = 1 : N
        x(j,ns) = alpha(j,1);
        y(j,ns) = beta(j,1);
        z(j,ns) = Gamma(j,1); 
    end
    
  end
  

