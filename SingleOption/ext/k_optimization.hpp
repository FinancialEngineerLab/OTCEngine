/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//	How to use
//	
//		1. Inherit this class 
//			with overriding definition of f 
//			with method parameters:	dx, max_itr, error_range
//		2. excute find_zero()
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

#ifndef	INCLUDE_FILE___kdb_optimization_hpp
#define	INCLUDE_FILE___kdb_optimization_hpp

//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include <kdb_BasicMath.hpp>
#include <kdb_LinearAlgebra.hpp>

struct	k_opt_NM_spec {		//	Newton Method spec
	double	x0;			//	initial seed
	double	dx;			//	delta x for derivativ calculation
	double	err_tol;	//	error tolerence
	int		max_itr;	//	if iteration goes beyond this number without finding zero, stop the routine and notify the failure
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class 	k_opt_NewtonMethod
{

public:

	k_opt_NewtonMethod();
	~k_opt_NewtonMethod();
	
	double	find_zero( k_opt_NM_spec spec );

protected:
	virtual	double f( double x );
};



/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// 	Program Name			: Levenberg Marquardt Method
// 	Version					: 0.0.0.0
// 	Auther					: Hangseob CHO
// 	Date					: 2009.01.09
// 	Modified by				: 
// 	Modified at				:
// 	Copyright				: KDB Quant team
// 	Description				: Levenberg Marquardt Method
// 	Related Doc. Name		: HS-081020-01 Reflective Barrier on LM.doc
/////////////////////////////////////////////////////////////////////////////////////////////////////*/

double	k_opt_reflexive_barrier( double a, double LB, double UB );
void	k_opt_reflexive_barrier_array( 
	
	//	input
	int dim_a, double *a, double *LB, double *UB,  
	
	//	output
	double *a_result
	
	);

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Levenberg-Marquardt Method 
//
//	Find such a that minimizes chi
//
//	y_original, y : R^dim_a -> R^nb_obs
//
//	y(a) = y_original( image of a in bounds )
//	
//
/////////////////////////////////////////////////////////////////////////////////////////////////////*/

class k_opt_LevenbergMarquardt
{
public:

	virtual void y_original( double *a, double *y_ele);
	
	void	find_Minimizing_Parameter ( 
				double 	*seed,				//	input:		initial seed
				double 	*a_result,			//	output:	minimizing parameter
				double 	*chi_min,			//	output:	minimized error
				int		b_print_process		//	input:		print lsqv
				);
			
	double 	chi( double *ay);
	void	y( double *a, double *y_ele);
	
	int 	dim_a;
	double	*LB;	
	double	*UB;
	//	LB[i]	<=	 a[i]	<=	UB[i]
	//	i = 0,	...	, dim_a -1
	
	int 	nb_obs;
	double	*y_obs;

};	

#endif