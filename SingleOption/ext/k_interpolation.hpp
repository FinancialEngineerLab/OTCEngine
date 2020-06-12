#ifndef KDB_INTERPOLATION_H
#define KDB_INTERPOLATION_H

#include <vector>
#include <iostream> //cout
using namespace std;
/**********************************************************************/
/* Program Name      : Interpolation library
/* version           : 0.0.0.0
/* author            : kicheon chang
/* date              : 2008.12.26
/* Modified by       : kicheon chang
/* Modified at       : 
/* Copyright         : KDB Quant team
/* Description       : Interpolation functions
/* Related Doc. Name : none
/**********************************************************************/
double k_intp_spline(vector<double>x, vector<double> y, double x_star, int iDim);

double k_intp_cubic_spline(vector<double> x, vector<double> y, double x_star);
double k_intp_linear_spline(vector<double> x, vector<double> y, double x_star);
double k_intp_polynomial(vector<double> x, vector<double> y, double x_star);
double k_intp_bi_linear(double x1, double x2, double y1, double y2, 
					  double z11, double z12, double z21, double z22, double x_star,double y_star);
double k_intp_linear(double x1, double x2, double y1, double y2, double x_star);

double k_intp_extrapolate(vector<double> x, vector<double> y, double x_star, int ext_method);

double	k_intp_bilin_interp( 
							   double xm, double ym,

							   double xd, double xu,
							   double yd, double yu,

							   //	V(x,y)'s
							   double Vdd, double Vdu,		
							   double Vud, double Vuu		
							   );

double	k_intp_bilin_interp(
							   double	xMin,	double	xMax, 	int 	nx,
							   double	yMin,	double	yMax, 	int 	ny,
							   double	**fxy,	double	x,	  	double 	y		);

double k_intp_linear( double x_min, double x_max, int nx, double *fx, double x_star );


#endif
