
#ifndef	INCLUDE_FILE__kdb_miscellaneous_hpp
#define	INCLUDE_FILE__kdb_miscellaneous_hpp

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
#include <ctime>

using namespace std;

double *** k_misc_new_double_array_HS( int I, int J, int K );

double	** k_misc_new_double_array_HS( int I, int J );

int		** k_misc_new_int_array_HS( int I, int J );

void 	k_misc_copy_array_HS( double *target, double *source, int I );

void 	k_misc_copy_array_HS( int *target, int *source, int I );

void 	k_misc_copy_array_HS( double **target, double **source, int I, int J );

void 	k_misc_swap_address_HS( double ****A, double ****B );

void 	k_misc_swap_address_HS( double ***A, double ***B );

void 	k_misc_swap_address_HS( double **A, double **B );

void	k_misc_delete_array_HS( double ***A, int I, int J );

void	k_misc_delete_array_HS( double **A, int I );

void	k_misc_delete_array_HS( int **A, int I );

void	k_misc_delete_array_HS( double *A );

void	k_misc_delete_array_HS( int *A );

int		k_misc_binsearch_lower_index( double *sorted_array, int num_ele, double key );

int		check_strict_increasing( int *a, int nb_a );

int		check_monotone_increasing( int *a, int nb_a );

double intp1d(double* tenor, double* curve, int nCol, double t);
double intp1d(double target,double *x,double *y,int min, int max);  //통합필요
double intp1d(double target, double* px, std::vector<double>& uold, int min, int max, int init_i);
double getforward(double* tenor, double* curve, double nCol, double t);
double getforward(const double tau, const double *r, const double *term, const long nb);  //통합필요

void ivlv(double spot, double **Lvol_kospi,
	double **Ivol_kospi,
	double* vol_term, long nb_vol_term, double* vol_strike, long nb_vol_strike,
	double* r_curve, double* r_term, long nb_r,
	double* div_curve, double* div_term, long nb_div);

void xlvolToIvol(double** Ivol,double* Ivolsurface_xl, long rows, long cols, long startindex);
string getFnameTimeStartingWith(string init_str);
#endif