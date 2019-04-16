/*
 *
 *  CauchyProdPower.h
 *
 *  This code determines the CauchyProduct of a vector a raised to a power P and stores it in b.
 *
 */

double cauchypower(double a[],double b[],int n,double P)
{
	int i=n, j; 
	double CauchyPow;
    
		CauchyPow = 0;
		for (j=0; j<=i-1; j++)
        {
			CauchyPow = (j*(P+1)+P-i)*a[j+1]*b[i-j] + CauchyPow;
		}
        CauchyPow = ((i+1)*P*a[i+1]*b[0] + CauchyPow)/((i+1)*a[0]);
    
	return CauchyPow;
}
