/*
 *
 *  CauchyProd.h
 *
 *  This code determines the CauchyProduct of two vectors a and b
 *
 */


double cauchyprod(double a[],double b[], int n)
{
	int i=n, j;
	double CauchyP;
    
		CauchyP = 0;
		for (j=0; j<=i; j++)
        {
			CauchyP = a[j]*b[i-j] + CauchyP;
		}

	return CauchyP;
}
