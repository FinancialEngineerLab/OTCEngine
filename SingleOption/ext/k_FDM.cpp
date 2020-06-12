
// << UPDATE >>
//   091030 if nx, ny < 3 the procedure stops

#include <kdb_FDM.hpp>

int	k_fdm_backward_induction_2F_FDM_OS_const_coeff(
			double **V_init, double **V_term, 
			double  dt,
			double  dx, int Nx,
			double  dy, int Ny,
			double a, double b, double c, double d, double e, double f)
/////////////////////////////////////////////////////////////////////////////////////
//
//	one time step backward induction
//	find V_init with given V_term with one time step from t_term to t_init
//
//	Underlying PDE:
//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
//		* a, b, c, d, e, f : constant
//
//	Gamma zero condition as boundary condition
//
//	Nx, Ny: numer of step I.E. V_init, V_term[Nx+1][Ny+1]
//
/////////////////////////////////////////////////////////////////////////////////////
{
	int	b_fail = 0;
	
	if(Nx<2 || Ny<2) return 1;

	int		i, j, I, J;
	double 	*A, *B, *C, *q;
	double	**V_btwn;		//	V between t_init and t_term ; V at t = ( t_init + t_term)/2
	
	V_btwn	= k_misc_new_double_array_HS( Nx+1, Ny+1);
	
	I = Nx - 1;
	J = Ny - 1;
	
	// first sweep  //////////////////////////////////////////////////////////////////////////////////////////////
	
	A = new double [I+1];
	B = new double [I+1];
	C = new double [I+1];
	q = new double [I+1];
	
	for ( j = 1; j <= J; j++ ) {
		
		for ( i = 1; i <= I; i++ ) {
			
			A[i] = a/(dx*dx) - b/(2.*dx);
			B[i] = -1./dt - 2.*a/(dx*dx) + c;
			C[i] = a/(dx*dx) + b/(2.*dx);
			q[i] = 	-0.5*f*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
					        /(4.*dx*dy) 
					-V_term[i][j]/dt;
		}
		// zero 2nd derivative on x_min and x_max
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		A[I] -=     C[I];	B[I] += 2.0*C[I];
		
		k_mx_tridiagonal_solver( A, B, C, q, 1, I);
		
		for ( i = 1; i <= I; i++ )
			V_btwn[i][j] = q[i];
		
		// zero 2nd derivative on x_min and x_max	
		V_btwn[0  ][j] = 2*V_btwn[1][j] - V_btwn[2  ][j];
		V_btwn[I+1][j] = 2*V_btwn[I][j] - V_btwn[I-1][j];
	}
	
	for ( i = 0; i <= I+1; i++ )  {
		V_btwn[i][0  ] = 2*V_btwn[i][1] - V_btwn[i][2  ];
		V_btwn[i][J+1] = 2*V_btwn[i][J] - V_btwn[i][J-1];
	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;
	
	// 2nd sweep ////////////////////////////////////////////////////////////////////////////////////////////////	
	
	A = new double [J+1];
	B = new double [J+1];
	C = new double [J+1];
	q = new double [J+1];
	
	for ( i = 1; i <= I; i++ ) {
		
		for ( j = 1; j <= J; j++ ) {
			
			A[j] = d/(dy*dy) - e/(2.*dy);
			B[j] = -1./dt - 2.*d/(dy*dy);
			C[j] = d/(dy*dy) + e/(2.*dy);
			q[j] = 	-0.5*f*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
					        /(4.*dx*dy) 
					-V_btwn[i][j]/dt;
		}
		// zero 2nd derivative on y_min and y_max
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		A[J] -=     C[J];	B[J] += 2.0*C[J];
		
		k_mx_tridiagonal_solver( A, B, C, q, 1, J);
		
		for ( j = 1; j <= J; j++ )
			V_init[i][j] = q[j];
		
		// zero 2nd derivative on y_min and y_max	
		V_init[i][0  ] = 2*V_init[i][1] - V_init[i][2  ];
		V_init[i][J+1] = 2*V_init[i][J] - V_init[i][J-1];
	}
	
	for ( j = 0; j <= J+1; j++ )  {
		V_init[0  ][j] = 2*V_init[1][j] - V_init[2  ][j];
		V_init[I+1][j] = 2*V_init[I][j] - V_init[I-1][j];
	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;	
	
	k_misc_delete_array_HS( V_btwn, Nx+1 );

	return b_fail;
}


int	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__Gamma_zero__FIRST_HALF(
	//	OUTPUT
	double **V_btwn, 
	//	INPUT
	double **V_term, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f)
	/////////////////////////////////////////////////////////////////////////////////////
	//
	//	FIRST HALF SWEEP of one time step backward induction
	//	find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Gamma zero condition as boundary condition
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int	b_fail = 0;
	
	if(Nx<2 || Ny<2) return 1;

	int		i, j, I, J;
	double 	*A, *B, *C, *q;

	I = Nx - 1;
	J = Ny - 1;

	// first sweep  //////////////////////////////////////////////////////////////////////////////////////////////

	A = new double [I+1];
	B = new double [I+1];
	C = new double [I+1];
	q = new double [I+1];

	for ( j = 1; j <= J; j++ ) {

		for ( i = 1; i <= I; i++ ) {

			A[i] = a/(dx*dx) - b/(2.*dx);
			B[i] = -1./dt - 2.*a/(dx*dx) + c;
			C[i] = a/(dx*dx) + b/(2.*dx);
			q[i] = 	-0.5*f*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
				/(4.*dx*dy) 
				-V_term[i][j]/dt;
		}
		// zero 2nd derivative on x_min and x_max
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		A[I] -=     C[I];	B[I] += 2.0*C[I];

		k_mx_tridiagonal_solver( A, B, C, q, 1, I);

		for ( i = 1; i <= I; i++ )
			V_btwn[i][j] = q[i];

		// zero 2nd derivative on x_min and x_max	
		V_btwn[0  ][j] = 2*V_btwn[1][j] - V_btwn[2  ][j];
		V_btwn[I+1][j] = 2*V_btwn[I][j] - V_btwn[I-1][j];
	}

	for ( i = 0; i <= I+1; i++ )  {
		V_btwn[i][0  ] = 2*V_btwn[i][1] - V_btwn[i][2  ];
		V_btwn[i][J+1] = 2*V_btwn[i][J] - V_btwn[i][J-1];
	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;

	return b_fail;
}

int	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__Gamma_zero__SECOND_HALF(
	//	OUTPUT
	double **V_init, 
	//	INPUT
	double **V_btwn, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f)
	/////////////////////////////////////////////////////////////////////////////////////
	//
	//	SECOND HALF SWEEP of one time step backward induction
	//	find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Gamma zero condition as boundary condition
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int	b_fail = 0;
	
	if(Nx<2 || Ny<2) return 1;

	int		i, j, I, J;
	double 	*A, *B, *C, *q;

	I = Nx - 1;
	J = Ny - 1;

	// 2nd sweep ////////////////////////////////////////////////////////////////////////////////////////////////	

	A = new double [J+1];
	B = new double [J+1];
	C = new double [J+1];
	q = new double [J+1];

	for ( i = 1; i <= I; i++ ) {

		for ( j = 1; j <= J; j++ ) {

			A[j] = d/(dy*dy) - e/(2.*dy);
			B[j] = -1./dt - 2.*d/(dy*dy);
			C[j] = d/(dy*dy) + e/(2.*dy);
			q[j] = 	-0.5*f*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
				/(4.*dx*dy) 
				-V_btwn[i][j]/dt;
		}
		// zero 2nd derivative on y_min and y_max
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		A[J] -=     C[J];	B[J] += 2.0*C[J];

		k_mx_tridiagonal_solver( A, B, C, q, 1, J);

		for ( j = 1; j <= J; j++ )
			V_init[i][j] = q[j];

		// zero 2nd derivative on y_min and y_max	
		V_init[i][0  ] = 2*V_init[i][1] - V_init[i][2  ];
		V_init[i][J+1] = 2*V_init[i][J] - V_init[i][J-1];
	}

	for ( j = 0; j <= J+1; j++ )  {
		V_init[0  ][j] = 2*V_init[1][j] - V_init[2  ][j];
		V_init[I+1][j] = 2*V_init[I][j] - V_init[I-1][j];
	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;

	return b_fail;
}

void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_lower_barrier_FIRST_HALF(
	//	OUTPUT
	double **V_btwn, 
	//	INPUT
	double **V_term, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f,
	double	*V_btwn_x_min, double *V_btwn_y_min)
	/////////////////////////////////////////////////////////////////////////////////////
	//
	//	FIRST HALF SWEEP of one time step backward induction
	//	find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Assume that if ( x <= Lower_barrier_of_x or y <= Lower_barrier_of_y ) 
	//		then V is replaced by given value V
	//
	//	Gamma zero condition on x_max and y_max
	//
	//	V_x_min : value along y with fixed x = x_min
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int		i, j, I, J;
	double 	*A, *B, *C, *q;
	I = Nx - 1;
	J = Ny - 1;

	// first sweep  //////////////////////////////////////////////////////////////////////////////////////////////

	A = new double [I+1];
	B = new double [I+1];
	C = new double [I+1];
	q = new double [I+1];

	for ( j = 1; j <= J; j++ ) {

		for ( i = 1; i <= I; i++ ) {

			A[i] = a/(dx*dx) - b/(2.*dx);
			B[i] = -1./dt - 2.*a/(dx*dx) + c;
			C[i] = a/(dx*dx) + b/(2.*dx);
			q[i] = 	-0.5*f*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
				/(4.*dx*dy) 
				-V_term[i][j]/dt;
		}
		//	zero 2nd derivative on x_max
		A[I] -=     C[I];	B[I] += 2.0*C[I];
		//	given value on x_min
		q[1] -=		A[1]*V_btwn_x_min[j];

		k_mx_tridiagonal_solver( A, B, C, q, 1, I);

		for ( i = 1; i <= I; i++ )
			V_btwn[i][j] = q[i];

		//	zero 2nd derivative on x_max
		V_btwn[I+1][j] = 2*V_btwn[I][j] - V_btwn[I-1][j];
		//	given value on x_min
		V_btwn[0][j] = V_btwn_x_min[j];
	}

	for ( i = 0; i <= I+1; i++ )  {

		//	zero 2nd derivative on y_max
		V_btwn[i][J+1] = 2*V_btwn[i][J] - V_btwn[i][J-1];
		//	given value on y_min
		V_btwn[i][0] = V_btwn_y_min[i];
	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;
}

void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_lower_barrier_SECOND_HALF(
	//	OUTPUT
	double **V_init, 
	//	INPUT
	double **V_btwn, 
	double  dt,			//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f,
	double	*V_init_x_min, double *V_init_y_min)
	/////////////////////////////////////////////////////////////////////////////////////
	//
	//	SECOND HALF SWEEP of one time step backward induction
	//	find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Assume that if ( x <= Lower_barrier_of_x or y <= Lower_barrier_of_y ) 
	//		then V is replaced by given value V
	//
	//	Gamma zero condition on x_max and y_max
	//
	//	V_x_min : value along y with fixed x = x_min
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int		i, j, I, J;
	double 	*A, *B, *C, *q;
	I = Nx - 1;
	J = Ny - 1;

	// 2nd sweep ////////////////////////////////////////////////////////////////////////////////////////////////	

	A = new double [J+1];
	B = new double [J+1];
	C = new double [J+1];
	q = new double [J+1];

	for ( i = 1; i <= I; i++ ) {

		for ( j = 1; j <= J; j++ ) {

			A[j] = d/(dy*dy) - e/(2.*dy);
			B[j] = -1./dt - 2.*d/(dy*dy);
			C[j] = d/(dy*dy) + e/(2.*dy);
			q[j] = 	-0.5*f*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
				/(4.*dx*dy) 
				-V_btwn[i][j]/dt;
		}
		//	zero 2nd derivative on y_max
		A[J] -=     C[J];	B[J] += 2.0*C[J];
		//	given value on y_min
		q[1] -=		A[1]*V_init_y_min[i];

		k_mx_tridiagonal_solver( A, B, C, q, 1, J);

		for ( j = 1; j <= J; j++ )
			V_init[i][j] = q[j];

		// zero 2nd derivative on y_max	
		V_init[i][J+1] = 2*V_init[i][J] - V_init[i][J-1];
		//	given value on y_min
		V_init[i][0] = V_init_y_min[i];
	}

	for ( j = 0; j <= J+1; j++ )  {

		//	zero 2nd derivative on x_max
		V_init[I+1][j] = 2*V_init[I][j] - V_init[I-1][j];
		//	given value on x_min
		V_init[0][j] = V_init_x_min[j];
	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;	
}


void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_upper_barrier_FISRT_HALF(
	//	OUTPUT
	double **V_btwn,
	//	INPUT
	double **V_term, 
	double  dt,		//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx, int Nx_UB,
	double  dy, int Ny,	int	Ny_UB,
	double a, double b, double c, double d, double e, double f,
	double  *V_btwn_x_UB, double *V_btwn_y_UB )
	/////////////////////////////////////////////////////////////////////////////////////
	//	
	//	FIRST HALF SWEEP of one time step backward induction
	//	find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Gamma zero condition as boundary condition except min_upper_barrier, where
	//		the values are given
	//
	//	Nx,    Ny    : number of asset steps (i.e. Nx+1, Ny +1 : number of asset points
	//	Nx_UB, Ny_UB : number of asset steps up to up barrier 
	//
	//	V_x_UB[j] : barrier value when x is at the barrier for	Ny_KO <= j <= Ny 
	//	
	//	V_y_UB[i] : barrier value when y is at the barrier	for Nx_KO <= i <= Ny 
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int		i, j, I, J, I_UB, J_UB;
	double 	*A, *B, *C, *q;

	I		= Nx	- 1;
	J		= Ny	- 1;
	I_UB	= Nx_UB - 1;
	J_UB	= Ny_UB - 1;

	// first sweep  //////////////////////////////////////////////////////////////////////////////////////////////

	A = new double [I+1];
	B = new double [I+1];
	C = new double [I+1];
	q = new double [I+1];

	//	step (1)
	for ( j = 1; j <= J_UB ; j++ ) 
	{
		for ( i = 1; i <= I; i++ ) 
		{
			A[i] = a/(dx*dx) - b/(2.*dx);
			B[i] = -1./dt - 2.*a/(dx*dx) + c;
			C[i] = a/(dx*dx) + b/(2.*dx);
			q[i] = 	-0.5*f*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
				/(4.*dx*dy) 
				-V_term[i][j]/dt;
		}
		// zero 2nd derivative on x_min and x_max
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		A[I] -=     C[I];	B[I] += 2.0*C[I];

		k_mx_tridiagonal_solver( A, B, C, q, 1, I);

		for ( i = 1; i <= I; i++ )
			V_btwn[i][j] = q[i];

		// zero 2nd derivative on x_min and x_max	
		V_btwn[0  ][j] = 2*V_btwn[1][j] - V_btwn[2  ][j];
		V_btwn[I+1][j] = 2*V_btwn[I][j] - V_btwn[I-1][j];
	}

	//	step (2)
	for ( j = J_UB + 1; j <= J; j++ ) 
	{
		for( i = 1; i <= I_UB; i++ )
		{
			A[i] = a/(dx*dx) - b/(2.*dx);
			B[i] = -1./dt - 2.*a/(dx*dx) + c;
			C[i] = a/(dx*dx) + b/(2.*dx);
			q[i] = 	-0.5*f*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
				/(4.*dx*dy) 
				-V_term[i][j]/dt;
		}
		//	zero 2nd derivative on x_min	
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		//	given value on x_UB
		q[I_UB] -= C[I_UB]*V_btwn_x_UB[j];

		k_mx_tridiagonal_solver( A, B, C, q, 1, I_UB);

		for ( i = 1; i <= I_UB; i++ )
			V_btwn[i][j] = q[i];

		//	zero 2nd derivative on x_min	
		V_btwn[0  ][j] = 2*V_btwn[1][j] - V_btwn[2  ][j];
		//	given value on x_UB
		V_btwn[I_UB+1][j] = V_btwn_x_UB[j];
	}

	//	step (3)
	for ( i = 0; i <= I+1; i++ )
		V_btwn[i][0] = 2*V_btwn[i][1] - V_btwn[i][2];

	//	step (4)
	for ( i = 0; i <= I_UB; i++ )  
		V_btwn[i][J+1] = 2*V_btwn[i][J] - V_btwn[i][J-1];
	
	//	step (5)
	V_btwn[I_UB + 1][J+1] =  V_btwn_x_UB[J+1];

	//	step (6)
	for ( i = I_UB+2; i <= I+1; i++ )
		V_btwn[i][J_UB+1] = V_btwn_y_UB[i];

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;
}

void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__min_upper_barrier_SECOND_HALF(
	//	OUTPUT
	double **V_init,
	//	INPUT
	double **V_btwn, 
	double  dt,		//	time step size of ( FIRST_HALF + SECOND_HALF ) ( not the half alone )
	double  dx, int Nx, int Nx_UB,
	double  dy, int Ny,	int	Ny_UB,
	double a, double b, double c, double d, double e, double f,
	double  *V_init_x_UB, double *V_init_y_UB )
	/////////////////////////////////////////////////////////////////////////////////////
	//	
	//	FIRST HALF SWEEP of one time step backward induction
	//	find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Gamma zero condition as boundary condition except min_upper_barrier, where
	//		the values are given
	//
	//	Nx,    Ny    : number of asset steps (i.e. Nx+1, Ny +1 : number of asset points
	//	Nx_KO, Ny_KO : number of asset steps up to KO barrier 
	//
	//	V_x_UB[j] : barrier value when x is at the barrier for	Ny_KO <= j <= Ny 
	//	
	//	V_y_UB[i] : barrier value when y is at the barrier	for Nx_KO <= i <= Ny 
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int		i, j, I, J, I_UB, J_UB;
	double 	*A, *B, *C, *q;

	I		= Nx	- 1;
	J		= Ny	- 1;
	I_UB	= Nx_UB - 1;
	J_UB	= Ny_UB - 1;

	// 2nd sweep ////////////////////////////////////////////////////////////////////////////////////////////////	

	A = new double [J+1];
	B = new double [J+1];
	C = new double [J+1];
	q = new double [J+1];

	//	step (1)
	for ( i = 1; i <= I_UB; i++ ) 
	{
		for ( j = 1; j <= J; j++ ) 
		{
			A[j] = d/(dy*dy) - e/(2.*dy);
			B[j] = -1./dt - 2.*d/(dy*dy);
			C[j] = d/(dy*dy) + e/(2.*dy);
			q[j] = 	-0.5*f*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
				/(4.*dx*dy) 
				-V_btwn[i][j]/dt;
		}
		// zero 2nd derivative on y_min and y_max
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		A[J] -=     C[J];	B[J] += 2.0*C[J];

		k_mx_tridiagonal_solver( A, B, C, q, 1, J);

		for ( j = 1; j <= J; j++ )
			V_init[i][j] = q[j];

		// zero 2nd derivative on y_min and y_max	
		V_init[i][0  ] = 2*V_init[i][1] - V_init[i][2  ];
		V_init[i][J+1] = 2*V_init[i][J] - V_init[i][J-1];
	}

	//	step (2)
	for ( i = I_UB + 1; i <= I; i++ )
	{
		for ( j = 1; j <= J_UB; j++ )
		{
			A[j] = d/(dy*dy) - e/(2.*dy);
			B[j] = -1./dt - 2.*d/(dy*dy);
			C[j] = d/(dy*dy) + e/(2.*dy);
			q[j] = 	-0.5*f*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
				/(4.*dx*dy) 
				-V_btwn[i][j]/dt;
		}
		//	zero 2nd derivative on y_min
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		//	given value on y_UB
		q[J_UB] -= C[J_UB]*V_init_y_UB[i];

		k_mx_tridiagonal_solver( A, B, C, q, 1, J_UB );

		for ( j = 1; j <= J_UB; j++ )
			V_init[i][j] = q[j];

		//	zero 2nd derivative on y_min
		V_init[i][0] = 2*V_init[i][1] - V_init[i][2];
		//	given value on y_UB
		V_init[i][J_UB+1] = V_init_y_UB[i];
	}

	//	step (3)
	for ( j = 0; j <= J+1; j++ )
		V_init[0  ][j] = 2*V_init[1][j] - V_init[2  ][j];

	//	step (4)
	for ( j = 0; j <= J_UB; j++ )
		V_init[I+1][j] = 2*V_init[I][j] - V_init[I-1][j];

	//	step (5)
	V_init[I+1][J_UB+1] = V_init_y_UB[I+1];

	//	step (6)
	for ( j = J_UB+2; j <= J+1; j++ )
		V_init[I_UB+1][j] = V_init_x_UB[j];

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;	
}

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
	double  *V_btwn_x_UB,  double *V_btwn_y_UB )
	/////////////////////////////////////////////////////////////////////////////////////
	//	
	//	IN PROGRESS
	//
	//	FIRST HALF SWEEP of one time step backward induction
	//		find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Boundary Condition:
	//		Gamma zero condition as boundary condition except min_upper_barrier and 
	//		min_lower_barrier, where the values are given
	//
	//	Nx,    Ny    : number of asset steps (i.e. Nx+1, Ny +1 : number of asset points
	//	Nx_UB, Ny_UB : number of asset steps up to up barrier 
	//
	//	V_x_UB[j] : barrier value when x is at the barrier for	Ny_KO <= j <= Ny 
	//	
	//	V_y_UB[i] : barrier value when y is at the barrier	for Nx_KO <= i <= Ny 
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int		i, j, I, J, I_UB, J_UB;
	double 	*A, *B, *C, *q;

	I		= Nx	- 1;
	J		= Ny	- 1;
	I_UB	= Nx_UB - 1;
	J_UB	= Ny_UB - 1;

	// first sweep  //////////////////////////////////////////////////////////////////////////////////////////////

	A = new double [I+1];
	B = new double [I+1];
	C = new double [I+1];
	q = new double [I+1];

	//	step (1)
	for ( j = 1; j <= J_UB ; j++ ) 
	{
		for ( i = 1; i <= I; i++ ) 
		{
			A[i] = a/(dx*dx) - b/(2.*dx);
			B[i] = -1./dt - 2.*a/(dx*dx) + c;
			C[i] = a/(dx*dx) + b/(2.*dx);
			q[i] = 	-0.5*f*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
				/(4.*dx*dy) 
				-V_term[i][j]/dt;
		}
		//  zero 2nd derivative on x_max
		A[I] -=     C[I];	B[I] += 2.0*C[I];
		//	given value on x_min
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		
		k_mx_tridiagonal_solver( A, B, C, q, 1, I);

		for ( i = 1; i <= I; i++ )
			V_btwn[i][j] = q[i];

		//	zero 2nd derivative on x_max
		V_btwn[I+1][j] = 2*V_btwn[I][j] - V_btwn[I-1][j];
		//	given value on x_min
		V_btwn[0][j] = V_btwn_x_min[j];
	}

	//	step (2)
	for ( j = J_UB + 1; j <= J; j++ ) 
	{
		for( i = 1; i <= I_UB; i++ )
		{
			A[i] = a/(dx*dx) - b/(2.*dx);
			B[i] = -1./dt - 2.*a/(dx*dx) + c;
			C[i] = a/(dx*dx) + b/(2.*dx);
			q[i] = 	-0.5*f*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
				/(4.*dx*dy) 
				-V_term[i][j]/dt;
		}
		//	given value on x_min	
		q[1] -=		A[1]*V_btwn_x_min[j];
		//	given value on x_UB
		q[I_UB] -= C[I_UB]*V_btwn_x_UB[j];

		k_mx_tridiagonal_solver( A, B, C, q, 1, I_UB);

		for ( i = 1; i <= I_UB; i++ )
			V_btwn[i][j] = q[i];

		//	given value on x_min
		V_btwn[0][j] = V_btwn_x_min[j];
		//	given value on x_UB
		V_btwn[I_UB+1][j] = V_btwn_x_UB[j];
	}

	//	step (3) : given value at y_min
	for ( i = 0; i <= I+1; i++ )
		V_btwn[i][0] = V_btwn_y_min[i];

	//	step (4) : Gamma zero at y_max
	for ( i = 0; i <= I_UB; i++ )  
		V_btwn[i][J+1] = 2*V_btwn[i][J] - V_btwn[i][J-1];

	//	step (5) : given value at y_max
	V_btwn[I_UB + 1][J+1] =  V_btwn_x_UB[J+1];

	//	step (6) : given value at y_UB
	for ( i = I_UB+2; i <= I+1; i++ )
		V_btwn[i][J_UB+1] = V_btwn_y_UB[i];

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;
}

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
	double  *V_init_x_UB,  double *V_init_y_UB )
	/////////////////////////////////////////////////////////////////////////////////////
	//	
	//	IN PROGRESS
	//
	//	FIRST HALF SWEEP of one time step backward induction
	//		find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Boundary Condition:
	//		Gamma zero condition as boundary condition except min_upper_barrier and 
	//		min_lower_barrier, where the values are given
	//
	//	Nx,    Ny    : number of asset steps (i.e. Nx+1, Ny +1 : number of asset points
	//	Nx_KO, Ny_KO : number of asset steps up to KO barrier 
	//
	//	V_x_UB[j] : barrier value when x is at the barrier for	Ny_KO <= j <= Ny 
	//	
	//	V_y_UB[i] : barrier value when y is at the barrier	for Nx_KO <= i <= Ny 
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int		i, j, I, J, I_UB, J_UB;
	double 	*A, *B, *C, *q;

	I		= Nx	- 1;
	J		= Ny	- 1;
	I_UB	= Nx_UB - 1;
	J_UB	= Ny_UB - 1;

	// 2nd sweep ////////////////////////////////////////////////////////////////////////////////////////////////	

	A = new double [J+1];
	B = new double [J+1];
	C = new double [J+1];
	q = new double [J+1];

	//	step (1)
	for ( i = 1; i <= I_UB; i++ ) 
	{
		for ( j = 1; j <= J; j++ ) 
		{
			A[j] = d/(dy*dy) - e/(2.*dy);
			B[j] = -1./dt - 2.*d/(dy*dy);
			C[j] = d/(dy*dy) + e/(2.*dy);
			q[j] = 	-0.5*f*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
				/(4.*dx*dy) 
				-V_btwn[i][j]/dt;
		}
		//	zero 2nd derivative on y_max
		A[J] -=     C[J];	B[J] += 2.0*C[J];
		//	given value on y_min
		q[1] -=		A[1]*V_init_y_min[i];

		k_mx_tridiagonal_solver( A, B, C, q, 1, J);

		for ( j = 1; j <= J; j++ )
			V_init[i][j] = q[j];

		// zero 2nd derivative on y_max	
		V_init[i][J+1] = 2*V_init[i][J] - V_init[i][J-1];
		//	given value on y_min
		V_init[i][0] = V_init_y_min[i];
	}

	//	step (2)
	for ( i = I_UB + 1; i <= I; i++ )
	{
		for ( j = 1; j <= J_UB; j++ )
		{
			A[j] = d/(dy*dy) - e/(2.*dy);
			B[j] = -1./dt - 2.*d/(dy*dy);
			C[j] = d/(dy*dy) + e/(2.*dy);
			q[j] = 	-0.5*f*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
				/(4.*dx*dy) 
				-V_btwn[i][j]/dt;
		}
		//	given value on y_min
		q[1] -=		A[1]*V_init_y_min[i];
		//	given value on y_UB
		q[J_UB] -= C[J_UB]*V_init_y_UB[i];

		k_mx_tridiagonal_solver( A, B, C, q, 1, J_UB );

		for ( j = 1; j <= J_UB; j++ )
			V_init[i][j] = q[j];

		//	given value on y_min
		V_init[i][0] = V_init_y_min[i];
		//	given value on y_UB
		V_init[i][J_UB+1] = V_init_y_UB[i];
	}

	//	step (3) : given value on x_min
	for ( j = 0; j <= J+1; j++ )
		V_init[0][j] = V_init_x_min[j];

	//	step (4) : Gamma zero on x_max
	for ( j = 0; j <= J_UB; j++ )
		V_init[I+1][j] = 2*V_init[I][j] - V_init[I-1][j];

	//	step (5) : given value on x_max
	V_init[I+1][J_UB+1] = V_init_y_UB[I+1];

	//	step (6) : given value on x_UB
	for ( j = J_UB+2; j <= J+1; j++ )
		V_init[I_UB+1][j] = V_init_x_UB[j];

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;	
}


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
	double	*V_btwn_x_bdry, double *V_btwn_y_bdry)
	/////////////////////////////////////////////////////////////////////////////////////
	//
	//	FIRST HALF SWEEP of one time step backward induction
	//	find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Assume that if ( x <= Lower_barrier_of_x or y <= Lower_barrier_of_y ) 
	//		then V is replaced by given value V
	//
	//	Gamma zero condition on x_max and y_max
	//
	//	V_x_min : value along y with fixed x = x_min
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int	b_fail = 0;
	
	if(Nx<2 || Ny<2) return 1;

	int		i, j, I, J, I0, J0, I1, J1;
	double 	*A, *B, *C, *q;
	I = Nx - 1;
	J = Ny - 1;

	I0 = (b_x_given_lower_bdry)? x_bdry_index: 0;
	J0 = (b_y_given_lower_bdry)? y_bdry_index: 0;
	I1 = I0+1;
	J1 = J0+1;

	// first sweep  //////////////////////////////////////////////////////////////////////////////////////////////

	A = new double [I+1];
	B = new double [I+1];
	C = new double [I+1];
	q = new double [I+1];

	for ( j = J1; j <= J; j++ ) {

		for ( i = I1; i <= I; i++ ) {

			A[i] = a/(dx*dx) - b/(2.*dx);
			B[i] = -1./dt - 2.*a/(dx*dx) + c;
			C[i] = a/(dx*dx) + b/(2.*dx);
			q[i] = 	-0.5*f*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
				/(4.*dx*dy) 
				-V_term[i][j]/dt;
		}
		//	zero 2nd derivative on x_max
		A[I] -=     C[I];	B[I] += 2.0*C[I];

		if(b_x_given_lower_bdry)
			q[I1] -=		A[I1]*V_btwn_x_bdry[j];     //  given value on x_min
		else
			B[1] += 2.0*A[1];	C[1] -=     A[1]; //  gamma zero on x_min

		k_mx_tridiagonal_solver( A, B, C, q, I1, I);

		for ( i = I1; i <= I; i++ )
			V_btwn[i][j] = q[i];

		//	zero 2nd derivative on x_max
		V_btwn[I+1][j] = 2*V_btwn[I][j] - V_btwn[I-1][j];

		if(b_x_given_lower_bdry)
			//	given value on x_min
			V_btwn[I0][j] = V_btwn_x_bdry[j];
		else 
			//	gamma zero on x_min
			V_btwn[0][j] = 2*V_btwn[1][j] - V_btwn[2][j];
	}

	for ( i = I0; i <= I+1; i++ )  {

		//	zero 2nd derivative on y_max
		V_btwn[i][J+1] = 2*V_btwn[i][J] - V_btwn[i][J-1];
		
		if(b_y_given_lower_bdry)
			//	given value on y_min
			V_btwn[i][J0] = V_btwn_y_bdry[i];
		else
			//  Gamma zero
			V_btwn[i][0  ] = 2*V_btwn[i][1] - V_btwn[i][2];		
	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;

	return b_fail;
}

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
	double	*V_init_x_bdry, double *V_init_y_bdry)
	/////////////////////////////////////////////////////////////////////////////////////
	//
	//	SECOND HALF SWEEP of one time step backward induction
	//	find V_btwn with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Assume that if ( x <= Lower_barrier_of_x or y <= Lower_barrier_of_y ) 
	//		then V is replaced by given value V
	//
	//	Gamma zero condition on x_max and y_max
	//
	//	V_x_min : value along y with fixed x = x_min
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int	b_fail = 0;
	
	if(Nx<2 || Ny<2) return 1;

	int		i, j, I, J, I0, J0, I1, J1;
	double 	*A, *B, *C, *q;
	I = Nx - 1;
	J = Ny - 1;

	I0 = (b_x_given_lower_bdry)? x_bdry_index: 0;
	J0 = (b_y_given_lower_bdry)? y_bdry_index: 0;
	I1 = I0+1;
	J1 = J0+1;

	// 2nd sweep ////////////////////////////////////////////////////////////////////////////////////////////////	

	A = new double [J+1];
	B = new double [J+1];
	C = new double [J+1];
	q = new double [J+1];

	for ( i = I1; i <= I; i++ ) {

		for ( j = J1; j <= J; j++ ) {

			A[j] = d/(dy*dy) - e/(2.*dy);
			B[j] = -1./dt - 2.*d/(dy*dy);
			C[j] = d/(dy*dy) + e/(2.*dy);
			q[j] = 	-0.5*f*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
				/(4.*dx*dy) 
				-V_btwn[i][j]/dt;
		}
		//	zero 2nd derivative on y_max
		A[J] -=     C[J];	B[J] += 2.0*C[J];
		//	given value on y_min
		if(b_y_given_lower_bdry)
			q[J1] -=		A[J1]*V_init_y_bdry[i];
		else
			B[1] += 2.0*A[1];	C[1] -=     A[1];

		k_mx_tridiagonal_solver( A, B, C, q, J1, J);

		for ( j = J1; j <= J; j++ )
			V_init[i][j] = q[j];

		// zero 2nd derivative on y_max	
		V_init[i][J+1] = 2*V_init[i][J] - V_init[i][J-1];
		if(b_y_given_lower_bdry)
			//	given value on y_min
			V_init[i][J0] = V_init_y_bdry[i];
		else
			V_init[i][0  ] = 2*V_init[i][1] - V_init[i][2  ];

	}

	for ( j = J0; j <= J+1; j++ )  {

		//	zero 2nd derivative on x_max
		V_init[I+1][j] = 2*V_init[I][j] - V_init[I-1][j];
		if(b_x_given_lower_bdry)
			//	given value on x_min
			V_init[I0][j] = V_init_x_bdry[j];
		else
			V_init[0  ][j] = 2*V_init[1][j] - V_init[2  ][j];

	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;	

	return b_fail;
}
	
	
	
/*
void	k_fdm__Bwd_Idctn__2F_FDM_OS_const_coeff__Gamma_zero(
	//	OUTPUT
	double **V_init,	
	double **V_btwn,
	//	INPUT
	double **V_term, 
	double  dt,
	double  dx, int Nx,
	double  dy, int Ny,
	double a, double b, double c, double d, double e, double f)
	/////////////////////////////////////////////////////////////////////////////////////
	//
	//	one time step backward induction
	//	find V_init with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + a*V_xx + f*V_xy + d*V_yy + b*V_x + e*V_y + c*V = 0
	//		* a, b, c, d, e, f : constant
	//
	//	Gamma zero condition as boundary condition
	//
	/////////////////////////////////////////////////////////////////////////////////////
{
	int		i, j, I, J;
	double 	*A, *B, *C, *q;

	I = Nx - 1;
	J = Ny - 1;

	// first sweep  //////////////////////////////////////////////////////////////////////////////////////////////

	A = new double [I+1];
	B = new double [I+1];
	C = new double [I+1];
	q = new double [I+1];

	for ( j = 1; j <= J; j++ ) {

		for ( i = 1; i <= I; i++ ) {

			A[i] = a/(dx*dx) - b/(2.*dx);
			B[i] = -1./dt - 2.*a/(dx*dx) + c;
			C[i] = a/(dx*dx) + b/(2.*dx);
			q[i] = 	-0.5*f*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
				/(4.*dx*dy) 
				-V_term[i][j]/dt;
		}
		// zero 2nd derivative on x_min and x_max
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		A[I] -=     C[I];	B[I] += 2.0*C[I];

		k_mx_tridiagonal_solver( A, B, C, q, 1, I);

		for ( i = 1; i <= I; i++ )
			V_btwn[i][j] = q[i];

		// zero 2nd derivative on x_min and x_max	
		V_btwn[0  ][j] = 2*V_btwn[1][j] - V_btwn[2  ][j];
		V_btwn[I+1][j] = 2*V_btwn[I][j] - V_btwn[I-1][j];
	}

	for ( i = 0; i <= I+1; i++ )  {
		V_btwn[i][0  ] = 2*V_btwn[i][1] - V_btwn[i][2  ];
		V_btwn[i][J+1] = 2*V_btwn[i][J] - V_btwn[i][J-1];
	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;

	// 2nd sweep ////////////////////////////////////////////////////////////////////////////////////////////////	

	A = new double [J+1];
	B = new double [J+1];
	C = new double [J+1];
	q = new double [J+1];

	for ( i = 1; i <= I; i++ ) {

		for ( j = 1; j <= J; j++ ) {

			A[j] = d/(dy*dy) - e/(2.*dy);
			B[j] = -1./dt - 2.*d/(dy*dy);
			C[j] = d/(dy*dy) + e/(2.*dy);
			q[j] = 	-0.5*f*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
				/(4.*dx*dy) 
				-V_btwn[i][j]/dt;
		}
		// zero 2nd derivative on y_min and y_max
		B[1] += 2.0*A[1];	C[1] -=     A[1];
		A[J] -=     C[J];	B[J] += 2.0*C[J];

		k_mx_tridiagonal_solver( A, B, C, q, 1, J);

		for ( j = 1; j <= J; j++ )
			V_init[i][j] = q[j];

		// zero 2nd derivative on y_min and y_max	
		V_init[i][0  ] = 2*V_init[i][1] - V_init[i][2  ];
		V_init[i][J+1] = 2*V_init[i][J] - V_init[i][J-1];
	}

	for ( j = 0; j <= J+1; j++ )  {
		V_init[0  ][j] = 2*V_init[1][j] - V_init[2  ][j];
		V_init[I+1][j] = 2*V_init[I][j] - V_init[I-1][j];
	}

	if (A)	delete[] A;
	if (B)	delete[] B;
	if (C)	delete[] C;
	if (q)	delete[] q;	
}




*/

