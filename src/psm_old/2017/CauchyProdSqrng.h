/*
 *
 *  CauchyProdSqrng.h
 *
 * This code determines the CauchyProduct for squaring a vector a
 *
 */


double cauchysquare(double a[], int n)
{
	int i=n, j;
	double CauchyPS; 
    
    CauchyPS = 0;
    if (i % 2 == 0)
    {
        for (j=0; j<=i/2-1; j++)
        {
			CauchyPS = a[j]*a[i-j] + CauchyPS;
		}
        CauchyPS = 2*CauchyPS + a[i/2]*a[i/2];
    }
    else
    {
        for (j=0; j<=(i-1)/2; j++)
        {
            CauchyPS = a[j]*a[i-j] + CauchyPS;
        }
        CauchyPS = 2*CauchyPS;
    }

	return CauchyPS;
}
