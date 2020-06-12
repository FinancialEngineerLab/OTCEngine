///////////////////////////////////////////////////////////////////////////////////
//                 KDB Probability  Header File
// 
//                 version 1.0.0 : 2009. 01. 08
// 
//
//                                   typed by ¤Ñ¤Ña
//                                        Jaemin Ahn
//
//
//////////////////////////////////////////////////////////////////////////////////

#ifndef KDB_Probability
#define KDB_Probability

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <limits>

#define SWITCH 3000	// when to switch to quadrature method


// PI
#ifndef PI
#define PI  3.14159265358979323846264338327950288419716939937510582
#endif


// 1/sqrt(2.0)
#ifndef ONE_OVER_SQRT_TWO 
#define ONE_OVER_SQRT_TWO     0.7071067811865475244008443621048490392848359376887
#endif


// 1/sqrt(2*PI)
#ifndef ONE_OVER_SQRT_TWO_PI
#define ONE_OVER_SQRT_TWO_PI  0.3989422804014326779399460599343818684758586311649
#endif


// sqrt(2*PI)
#ifndef SQRT_TWO_PI
#define SQRT_TWO_PI           2.506628274631000502415765284811045253006986740610
#endif



static const double cof[28] = {-1.3026537197817094, 6.4196979235649026e-1,
						 1.9476473204185836e-2, -9.561514786808631e-3, -9.46595344482036e-4,
						 3.66839497852761e-4, 4.2523324806907e-5, -2.02785781122534e-5,
						 -1.624290004647e-6, 1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8,
						 6.529054439e-9, 5.059343495e-9, -9.91364156e-10, -2.27365122e-10,
						 9.6467911e-11, 2.394038e-12, -6.886027e-12, 8.94487e-13, 3.13092e-13,
						 -1.12708e-13, 3.81e-16, 7.106e-15,-1.523e-15, -9.4e-17, 1.21e-16, -2.8e-17};
static const int ncof = 28;


double k_pr_cdf_T_two_tailed(double t, double dof);
//the two-tailed cdf of student's t-distribution, prob(<|t|), t>0

double k_pr_cdf_F(double f, double dof1, double dof2);
// cumulative distribution function of F distribution, prob(0<F<f)

double k_pr_ln_gamma(const double xx);
// ln[ Gamma(xx) ] for xx>0


double k_pr_CP(double x, int k);
// Cumulative Poisson Probability Function P_x(<k) integer k >=1
// x : expected mean number 

double k_pr_chi_P(double chis, int dof);
// the probability that the observed chi-square for a correct model
//  should be less than a value chis

double k_pr_chi_Q(double chis, int dof);
// the probability that the observed chi-square will exceed the value 
//   chis by chance even for a correct model.

double k_pr_Npdf(double mu, double sigma, double x);
// Probability density function for N(mu,sig)

double k_pr_CN(double x);  // mu = 0, sig = 1
//  Cummulative Univariate Normal Function  
//  J. Hart(1968) Computer Approximations, Wiley. Algorithm 5666 for the error function
//  modified by G. West(2004)
//  : accurate to double precision throughout the real line

double k_pr_CBN(double x, double y, double rho); // mu = 0, sig = 1
// computing bivariate cumulative normal probabilities                        
// A. Genz (2004) 'Numerical computation of rectangular bivariate and trivariate normal and
//                 t probabilities', Statistics and Computing 14 251-260
// accurate to double precision

double k_pr_ICN(double p);   //  0 < u < 1
// The Inverse cumulative normal distribution function
// by Peter J. Acklam, University of Oslo, Statistics Division.
//
// URL: http://www.math.uio.no/~jacklam/notes/invnorm





/* minor functions for evaluating t-distribution, f-distribution */
double k_pr_ibeta(double a, double b, double x);
// incomplete beta function I_x(a,b) for positive a and b, and x between 0 and 1

double k_pr_cf_beta(const double a, const double b, const double x);
// continued fraction for incomplete beta function by modified Lentz's method.

double k_pr_quad_ibeta(double a, double b, double x);
// Incomplete beta by quadrature. Returns I_x(a,b)




/* minor functions for evaluationg Possion distribution, Chi-square distribution */
void k_pr_gser(double *gamser, double a, double x, double *gln);
// Returns the incomplete gamma function P(a,x) evaluated by its series representation 
// as "gamser".  Also returns ln(Gamma(a)) as "gln"
//   x >= 0

void k_pr_gcf(double *gammcf, double a, double x, double *gln);
// Returns the incomplete gamma function Q(a,x) evaluated 
// by its continued fraction representation as "gamser".
// Also returns ln(Gamma(a)) as "gln"

double k_pr_igammap(double a, double x);
// Returns the incomplete gamma function P(a,x) 
// x >= 0, a > 0 

double k_pr_igammaq(double a, double x);
// Returns the incomplete gamma function Q(a,x)=1-P(a,x) 
// x >= 0, a > 0 




double k_pr_erfccheb(double z);
// erf using Chebyshev coefficients
// z >= 0

double k_pr_erf(double x);
// erf(x) error function

double k_pr_erfc(double x);
// erfc(x) complementary error function

double k_pr_inverfc(double p);  
// inverse of complementary error function
// 0 < p < 2

double k_pr_inverf(double p);
// inverse of error function
// -1 < p < 1

double k_pr_LNpdf(double mu, double sig, double x);
// probability density function of Log-normal distribution
// x >= 0


double k_pr_LNcdf(double mu, double sig, double x);
// cumulative distribution function of Log-normal distribution
// x >= 0

double k_pr_LNinvcdf(double mu, double sig, double p);
// inverse cumulative distribution function of Log-normal distribution
// 0 < p < 1


#endif





