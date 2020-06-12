#include "k_BasicMath.hpp"



/************************************************************************/
//	Return the nearest whole number
//		roundNearestWhole(2.4) => 2.0
//		roundNearestWhole(2.5) => 3.0
//		roundNearestWhole(2.6) => 3.0
//		roundNearestWhole(-2.4) => -2.0
//		roundNearestWhole(-2.5) => -2.0
//		roundNearestWhole(-2.6) => -3.0
/************************************************************************/
double k_bmath_roundNearestWhole(const double number)
{

	double fractpart, intpart;
	
	fractpart = modf(number , &intpart);

    return fabs(fractpart)>=0.5? intpart+k_bmath_sgn(number) : intpart ;
}


/************************************************************************/
// Returns the smallest integer value that is not less than a.
// Returns the same number with ceil() in <math.h>, but in int type
/************************************************************************/
int k_bmath_roundDown( double a)
{
	if ( a>=0 ) return (int)a;
	else if ( -a -(int)(-a) ==0) return -(int)(-a);
	else return -(int)(-a) -1;
}


/************************************************************************/
// Returns the largest integral value that is not greater than a.
// Returns the same number with floor() in <math.h>, but in int type
/************************************************************************/
int k_bmath_roundUp( double a)
{
	return -k_bmath_roundDown( -a);
}

/************************************************************************/
// Extracts the sign of a real number
/************************************************************************/
int k_bmath_sgn( double a)
{	
	if (a >=0.0) return 1;
	else return -1;	
}


int k_bmath_sgn( int a)
{	
	if (a >=0) return 1;
	else return -1;
}


double k_bmath_bound( double LB, double x, double UB )
{
	return min( max( LB, x), UB);
}

int	k_bmath_into_proper_integer( double A )
//============================================================================================
//	restore integer value stored as double type into value as integer type
//
//	the danger of (int)A is
//	
//	int		a_int		= 12;
//	double	a_double	= (double)a_int;
//
//	a_double may be either 12.00000000000001 or 11.9999999999999
//	the danger of returning the value integer by (int)a_double is
//	It may return 11 instead of 12.
//
//	So when changing type of value that is integer into int, we must use the below function
//		*added by Hangseob at 091105
//============================================================================================
{
	double	margin = 0.0000001;

	int	B = (int) A;
	double C = (double) B;

	if      ( A - C >  0.5 )  return B + 1;
	else if ( C - A >  0.5 )  return B - 1; 
	else				      return B;

	//	B can be either 1) the intended integer value or 2) smaller by 1.
	//	if 1), A and C must very close
	//	if 2), C must be smaller that A by about 1
}



double	double_max_in_array_HS ( double *a, int n, int m )
{
	int			i;
	double		maximum;

	maximum	= a[n];

	for ( i = n+1; i <= m; i++ )
		if	(a[i] > maximum)	maximum = a[i];

	return maximum;
}

double	double_min_in_array_HS ( double *a, int n, int m )
{
	int			i;
	double		minimum;

	minimum	= a[n];

	for ( i = n+1; i <= m; i++ )
		if	(a[i] < minimum)	minimum = a[i];

	return minimum;
}


int	power_HS( int m, int n )
{
	int result = 1;

	int	i;

	for ( i = 1; i <= n; i++ )
		result *= m;

	return 	result;
}

double k_power(double a, double x)
  // return a^x (a > 0.0);
{
      if (a < 0.0)  return(-1000000000000000.0); // error value
      if (a == 0.0)	return(0.0);

	  return( exp(x*log(a)) );

}




