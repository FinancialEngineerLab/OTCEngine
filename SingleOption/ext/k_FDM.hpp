
#ifndef	INCLUDE_FILE__kdb_FDM_hpp
#define	INCLUDE_FILE__kdb_FDM_hpp

#include <kdb_miscellaneous.hpp>
#include <kdb_LinearAlgebra.hpp>

int	k_fdm_backward_induction_2F_FDM_OS_const_coeff(
			double **V_init, double **V_term, 
			double  dt,
			double  dx, int Nx,
			double  dy, int Ny,
			double a, double b, double c, double d, double e, double f);
//	one time step backward induction
//	find V_init with given V_term with one time step from t_term to t_init			
//	Gamma zero condition as boundary condition

int k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__Gamma_zero__FIRST_HALF(
	//	OUTPUT
	double **V_btwn, 
	//	INPUT
	double **V_term, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f);

int	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__Gamma_zero__SECOND_HALF(
	//	OUTPUT
	double **V_init, 
	//	INPUT
	double **V_btwn, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f);

void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_lower_barrier_FIRST_HALF(
	//	OUTPUT
	double **V_btwn, 
	//	INPUT
	double **V_term, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f,
	double	*V_btwn_x_min, double *V_btwn_y_min);

void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_lower_barrier_SECOND_HALF(
	//	OUTPUT
	double **V_init, 
	//	INPUT
	double **V_btwn, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f,
	double	*V_init_x_min, double *V_init_y_min);

void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_upper_barrier_FISRT_HALF(
	//	OUTPUT
	double **V_btwn,
	//	INPUT
	double **V_term, 
	double  dt,		//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx, int Nx_UB,
	double  dy, int Ny,	int	Ny_UB,
	double a, double b, double c, double d, double e, double f,
	double  *V_btwn_x_UB, double *V_btwn_y_UB );

void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_upper_barrier_SECOND_HALF(
	//	OUTPUT
	double **V_init,
	//	INPUT
	double **V_btwn, 
	double  dt,		//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx, int Nx_UB,
	double  dy, int Ny,	int	Ny_UB,
	double a, double b, double c, double d, double e, double f,
	double  *V_init_x_UB, double *V_init_y_UB );

void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_upper_and_lower_barrier_FISRT_HALF(
	//	OUTPUT
	double **V_btwn,
	//	INPUT
	double **V_term, 
	double  dt,		//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx, int Nx_UB,
	double  dy, int Ny,	int	Ny_UB,
	double a, double b, double c, double d, double e, double f,
	double	*V_btwn_x_min, double *V_btwn_y_min,
	double  *V_btwn_x_UB,  double *V_btwn_y_UB );

void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_upper_and_lower_barrier_SECOND_HALF(
	//	OUTPUT
	double **V_init,
	//	INPUT
	double **V_btwn, 
	double  dt,		//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx, int Nx_UB,
	double  dy, int Ny,	int	Ny_UB,
	double a, double b, double c, double d, double e, double f,
	double	*V_init_x_min, double *V_init_y_min,
	double  *V_init_x_UB,  double *V_init_y_UB );

int	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_lower_barrier_FIRST_HALF(
	//	OUTPUT
	double **V_btwn, 
	//	INPUT
	double **V_term, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f,
	int		b_x_given_lower_bdry, int b_y_given_lower_bdry,
	int		x_bdry_index, int y_bdry_index,
	double	*V_btwn_x_min, double *V_btwn_y_min);

int	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_lower_barrier_SECOND_HALF(
	//	OUTPUT
	double **V_init, 
	//	INPUT
	double **V_btwn, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f,
	int		b_x_given_lower_bdry, int b_y_given_lower_bdry,
	int		x_bdry_index, int y_bdry_index,
	double	*V_init_x_min, double *V_init_y_min);

#endif