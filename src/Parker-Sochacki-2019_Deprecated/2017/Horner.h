/*
 *
 * Horner.h
 *
 *  This determines the value of a polynomial at a given t using Horner's algorithm. 
 * 
 */

double horner(double c[], double t, int m)
{
	int i,n;
    double HornerSum;
    
    n = m;
    HornerSum = c[n-1] + t*c[n];
	for ( i = n-2; i >= 0; i--)
    {
		HornerSum = c[i] + t*HornerSum;
	}
    
	return HornerSum;
}
