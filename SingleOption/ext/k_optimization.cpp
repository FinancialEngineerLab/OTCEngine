///////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//	How to use
//	
//		1. Inherit this class 
//			with overriding definition of f 
//			with method parameters:	dx, max_itr, error_range
//		2. excute find_zero()
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <kdb_optimization.hpp> 

k_opt_NewtonMethod::k_opt_NewtonMethod()
{
}

k_opt_NewtonMethod::~k_opt_NewtonMethod()
{
}

double	k_opt_NewtonMethod::f( double x )
{
	printf("error!!! dummy function has been called");
	
	exit(0);
	
	return 0.;
}

double	k_opt_NewtonMethod::find_zero( k_opt_NM_spec spec )
{

//	find x such that 	f(x) = 0
//	by Newton method with initial x = x0
//	until abs( f(  x found )  ) <= error_range
//	If iteration goes beyond max_itr, stop the routine and notify the failure

	double 	x_old, x_new, error; 
	double 	deriv_f;
	
	int		nb_itr;		
	
	x_old 	= spec.x0;
	error 	= spec.err_tol + 1.0;	//	so that the below while loop will do the first routine
	
	nb_itr =0;
	
	while ( error > spec.err_tol && nb_itr <= spec.max_itr) {
		
		double	f_x_old = f(x_old);
		deriv_f = ( f(x_old + spec.dx) - f_x_old )/spec.dx;
		x_new = x_old - f_x_old/deriv_f;
		
		error = fabs( f(x_new) );
		
		nb_itr++;
		x_old = x_new;
	}
	
	if ( nb_itr == spec.max_itr + 1 ) {
		printf("\n Newton Method Error: can't find the solution \n");
		exit(1);
	}
	
	return x_new;
}


/*//////////////////////////////////////////////////////////////////////////////////////
//
//	Find mirror image of a that is inside [LB,UB] (LB: Lower Bound, UB: Upper Bound )
//
/*//////////////////////////////////////////////////////////////////////////////////////

double	k_opt_reflexive_barrier( double a, double LB, double UB )
{
	if 		( a > UB )	{
		int 	multi;
		double 	reminder;
		multi 		= (int)floor( (a - UB)/(UB - LB) );
		reminder 	=  a - ( UB + (double)multi * (UB - LB) );
		if ( multi%2 == 0 )	
			return UB - reminder;
		else
			return LB + reminder;
	}
	else if ( a < LB ) {
		return k_opt_reflexive_barrier( 2.0*LB - a, LB, UB );
	}
	else {
		return a;
	}	
}


/*//////////////////////////////////////////////////////////////////////////////////////
//
//	Find mirror image of (a[0], a[1], ... , a[dim_a-1] )
//	that is inside [LB[0],UB[0]]*[LB[1],UB[1]]* ... *[LB[dim_a-1],UB[dim_a-1]] 
//
//	dim_a : dimension of a
/*//////////////////////////////////////////////////////////////////////////////////////

void	k_opt_reflexive_barrier_array( 
	
	//	input
	int dim_a, double *a, double *LB, double *UB,  
	
	//	output
	double *a_result
	
	)
{
	int i;
	
	for ( i = 0; i < dim_a; i++ ) {
		a_result[i] = k_opt_reflexive_barrier( a[i], LB[i], UB[i] );		
	}

}

//	Levenberg-Marquardt Method 
//	y is discrete function
//	use reflective barrier to find minimizing a[i]'s within each bounds

void k_opt_LevenbergMarquardt::y( double *a, double *y_ele)
{
	int i;
	
	double *a_IB;	//	a in bounds
	a_IB = new double [dim_a];
	for ( i = 0; i < dim_a; i++)
		a_IB[i] = k_opt_reflexive_barrier( a[i], LB[i], UB[i] );
	
	y_original( a_IB, y_ele );
	
	if ( a_IB ) delete[] a_IB;
}


void k_opt_LevenbergMarquardt::find_Minimizing_Parameter(
			double 	*seed,		//	input:		initial seed
			double 	*a_result,				//	output:	minimizing parameter
			double 	*chi_min,			//	output:	minimized error
			int		b_print_process		//	input:		print lsqv
			)
//==============================================================================
//
//	* If k_mx_Alpha_System_Solver fails, set chi_min = -100.
//
//==============================================================================
{
	
	//	This function is implemented by DR. Ahn
		
	int    i,j,k;

	double d_success;
	
	//{	declare variables and memory allocation

	double	*a_init;
	a_init 			= (double*	)malloc(dim_a*sizeof(double		));
	for(i= 0;i< dim_a; i++)	
		a_init[i] 	= seed[i];
	
	double	**alpha;
	alpha 			= (double** )malloc(dim_a*sizeof(double*	));
	for (i=0; i< dim_a; i++)
		alpha[i] 	= (double* 	)malloc(dim_a*sizeof(double		));
	
	double	*beta;
	beta			= (double* 	)malloc(dim_a*sizeof(double		));
	
	double	*delta;
	delta 			= (double* 	)malloc(dim_a*sizeof(double		));
	
	int    count;
	double 	lambda;
	double 	lsqv1, lsqv2;	//	for storing and comparing least-square values
	
	double 	**ey;										 
	ey 				= (double**	)malloc(dim_a*sizeof(double*	));		
	for (i=0; i< dim_a; i++)								
		ey[i] 		= (double*	)malloc(nb_obs*sizeof(double	));	
	
	double	*ay;
	ay 				= (double*	)malloc(nb_obs*sizeof(double	));
	
	double 	*a_next;
	a_next 			= (double*	)malloc(dim_a*sizeof(double		));
	
	/////////////////////////////////////////////////////////////////////////////////////////
	//	modify: 2008-03-26
	//		check the lsqv decrease in the lastest 5 
	//			if the sum of deacrese is smaller than certain number stop the process
	
	
	
	//}
	
	//	set the method parameters
	double eps = 0.0000001;
	double TOL = 0.0001;

	//   	store y(x,a)
	y(a_init,ay);          
	   
	// 	store epsilon-varied y(x,  a+epislion e_j)   
	for(j=0; j < dim_a; j++) {
		
		a_init[j] +=  eps;
		
		y(a_init,ey[j]);	

		a_init[j] -= eps;
	}

	lambda 	= 0.001;                         // scaling factor 
    lsqv1 	= sqrt(chi(ay)/(double)nb_obs);                           

	if	(b_print_process) {
		// printf("\n");
		
		// for ( i = 0; i < dim_a; i++ )
			// printf 
		
		printf("\n first lsqv = %.15f \n", lsqv1);
		
	}
	
	lsqv2 	= 0.0;
    
	//  making  beta
	for(i=0; i < dim_a; i++){
		beta[i] = 0.0;
		for(j=0; j < nb_obs; j++){
			beta[i] += (y_obs[j]-ay[j])*(ey[i][j]-ay[j])/eps/y_obs[j]/y_obs[j];
		}
	}
	
	// making alpha
	for(i=0; i < dim_a; i++){
		for(j=i; j < dim_a; j++){
			alpha[i][j] = 0.0;
			for(k=0; k < nb_obs; k++){
				alpha[i][j] += (ey[i][k]-ay[k])*(ey[j][k]-ay[k])
								/eps/eps/y_obs[k]/y_obs[k];
			}
			alpha[j][i] = alpha[i][j];
	   }
    }
	
	// diagonal setting
	for(i=0; i < dim_a; i++)
		alpha[i][i] *= (1.0+lambda);

    count = 0;

	for(; (sqrt(lsqv1/nb_obs) >= TOL) && (count < 20) ; )	{
	  
		d_success = k_mx_Alpha_System_Solver(dim_a, alpha,delta,beta);
		
		if(d_success < 0.5) break;

		for(i=0; i < dim_a; i++) {
			a_next[i] = a_init[i] + delta[i];	  
		}
		
        y(a_next,ay);
		
		lsqv2 = sqrt(chi(ay)/(double)nb_obs);
	   
		if (lsqv2 >= lsqv1){
	         
			for (i=0; i < dim_a; i++){
				alpha[i][i] *= (1.0+10.0*lambda)/(1.0+lambda);
			}
			lambda *= 10.0;
			count++;
		}
		else
		{
			
			if	(b_print_process) { //	 put it into 
				
				printf("\n %d   ",count);
				
				for(i=0; i < dim_a; i++) 
					printf("%.4f   ", a_next[i]);
				
				printf("%.10f   ", lsqv1-lsqv2);
				printf("lsqv %.10f \n",lsqv2);
			}	
			
			////////////////////////////////////////////////////////////////////////////////////////////
			//	hangseob's stop algorithm
			if ( lsqv1 - lsqv2 <= 1./1000./1000.) break;
			////////////////////////////////////////////////////////////////////////////////////////////
			
			count =0;
		   
			lsqv1 = lsqv2;
			
			for(i=0; i < dim_a; i++)
				a_init[i] = a_next[i];
		   
			y(a_init,ay);
	        
			for(j=0; j < dim_a; j++) {
	            a_init[j] +=  eps;
                y(a_init,ey[j]);		  
		        a_init[j] -= eps;
	        }
			
			lambda /= 10.0;
		   
			for(i=0; i < dim_a; i++){
				beta[i] = 0.0;
				
				for(j=0; j < nb_obs; j++){
					beta[i]	+= 	(y_obs[j]-ay[j])*(ey[i][j]-ay[j])
								/eps/y_obs[j]/y_obs[j];
				}
			}
			
			for(i=0; i < dim_a; i++){
				for(j=i; j < dim_a; j++){
					
					alpha[i][j] = 0.0;
					
					for(k=0; k < nb_obs; k++){
						alpha[i][j]	+= 	(ey[i][k]-ay[k])*(ey[j][k]-ay[k])
										/eps/eps/y_obs[k]/y_obs[k];
					}
					
					alpha[j][i] = alpha[i][j];
				}
            }
			
			for(i=0; i < dim_a; i++)
				alpha[i][i] *= (1.0+lambda);
			
	    }
    }
	
	for(i=0; i < dim_a; i++) 
		a_result[i] = k_opt_reflexive_barrier( a_next[i], LB[i], UB[i] );
   
	//	return resulting    sqrt( chi)
	y(a_result,ay);
	*chi_min	= chi(ay);	
	
	if(d_success < 0.5) *chi_min = -100.;

	//	free	dynamically allocated memories
	if	(a_init)	free( a_init);
	for	(i=0; i < dim_a; i++) 
		if	(alpha[i]) 	free(alpha[i]);
	if	(alpha) 	free(alpha);	
	if	(beta) 		free(beta);
	if	(delta) 	free(delta);
	for	(i=0; i< dim_a; i++)
		if	(ey[i]) 	free(ey[i]);
	if	(ey) 		free(ey);	
	if	(ay) 		free(ay);
	if	(a_next) 	free(a_next);
		
}
 
double k_opt_LevenbergMarquardt::chi( double *ay)
{
	int 	i;
	double 	temp = 0.0;
	double 	error;
	
	for(i=0; i< nb_obs; i++) {
	
		error 	=  ( y_obs[i] - ay[i] ) / y_obs[i];
		temp 	+= error*error;
	
	}
	
   return temp;

} 

void k_opt_LevenbergMarquardt::y_original( double *a, double *y_ele)
{
	printf("error!!! dummy function has been called");
	exit(0);
}