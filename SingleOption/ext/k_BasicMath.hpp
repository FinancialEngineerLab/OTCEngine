#ifndef KDB_BASICMATH_H
#define KDB_BASICMATH_H

#include <math.h>
#include <algorithm>

using namespace std;

/**********************************************************************/
/* Program Name      : BasicMath library
/* version           : 0.1.0.1
/* author            : Jihyun Lee
/* date              : 2009. 1. 7
/* Modified by       : Jaemin Ahn
/* Modified at       : 2011. 1. 18
/* Copyright         : KDB Quant team
/* Description       : round, sort, etc
/* Related Doc. Name : none
/**********************************************************************/

double k_bmath_roundNearestWhole(double number);
int k_bmath_roundDown( double a);
int k_bmath_roundUp( double a);
int k_bmath_sgn( double a);
int k_bmath_sgn( int a);
double k_bmath_bound( double LB, double x, double UB );
int	k_bmath_into_proper_integer( double A );

double	double_max_in_array_HS ( double *a, int n, int m );
double	double_min_in_array_HS ( double *a, int n, int m );
int	power_HS( int m, int n );

double k_power(double a, double x);
// return a^x (a > 0.0);

#endif