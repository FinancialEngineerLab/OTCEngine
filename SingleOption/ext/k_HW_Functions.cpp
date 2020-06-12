#include "kdb_HW_Functions.hpp"

double k_phi1( 
			double T,           // time
			double dt,          // time-step for finite-difference
			double *para,       // HW parameters para[0]: mean reversion, para[1]: volatility
			k_cZeroCurve *ZC    // Zero Curve Structure
			)
// HW 1 Factor short rate mean level function r(T)t = x(T) + phi1(T)
{
	double A;

	double y, y_dt;
	double fM;

	y 	 =  log( (*ZC).DF(T) );
	y_dt =  log( (*ZC).DF(T + dt ) );

	fM =  - ( y_dt - y )/dt;

	A = 1.0 - exp(-para[0]*T);

	return fM + para[1]*para[1]/(2.0*para[0]*para[0])*A*A;

}

double k_PtT1(
			double t,         // start time t
			double T,         // maturity
			double dt,        // time-step for finite-difference
			double r,         // short rate at t
			double *para,     // HW paraemters para[0]: mean reversion, para[1]: volatility
			k_cZeroCurve *ZC    // Zero Cureve Structure
			)
// zero coupon bond price at time t maturing at T with unit face value 1 
//  under HW 1 factor model
{
	double A,B;

	double y, y_dt;
	double fM;

	y 	 =  log( (*ZC).DF(t) );
	y_dt =  log( (*ZC).DF(t + dt ) );

	fM =  - ( y_dt - y )/dt;

	B = (1.0 - exp(-para[0]*(T-t)))/para[0];

	A = (*ZC).DF(T)/(*ZC).DF(t) * exp(B*fM-para[1]*para[1]/para[0]/4.0*(1.0-exp(-2*para[0]*t))*B*B );


	return A*exp(-B*r);

}


double	k_Vddt1( 
			  double ddt,   // T-t
			  double *para  // HW parameters para[0]: mean reversion, para[1]: volatility 
			  )
			  // stochastic variable: R(t,T) = integral of r(u) between t and T
			  // return V(t,T) = variance of R(t,T) conditional of F_t under HW 1 Factor model
{
	double ea;
	double sum;

	ea 		= exp(-para[0]*ddt);

	sum	= ddt +2.0/para[0]*ea - 0.5/para[0]*ea*ea -1.5/para[0];

	return	  para[1]*para[1]/para[0]/para[0]*sum;
}

double	k_phi2( 
			 double T,           // time
			 double dt,          // time-step for finite-difference
			 double *para,       // HW parameters para[0]: mean reversion 1, para[1]: mean reversion 2
			                     //               para[2]: volatility 1,     para[3]: volatility 2
								 //               para[4]: correlation
			 k_cZeroCurve *ZC    // Zero Curve Structure
			 )
// HW 2 Factor short rate mean level function r(T) = x(T) + y(T) + phi2(T)
{
	double A, B;

	double y, y_dt;
	double fM;

	y 	 =  log( (*ZC).DF(T) );
	y_dt =  log( (*ZC).DF(T + dt ) );

	fM =  - ( y_dt - y )/dt;

	A = 1.0 - exp(-para[0]*T);
	B = 1.0 - exp(-para[1]*T);

	return fM + para[2]*para[2]/(2.0*para[0]*para[0])*A*A
		+ para[3]*para[3]/(2.0*para[1]*para[1])*B*B
		+ para[4]*(para[2]*para[3])/(para[0]*para[1])*A*B;
}

double	k_Vddt2(
			  double ddt,   // T-t
			  double *para  // HW paraemters para[0]: mean reversion 1, para[1]: mean reversion 2
							//               para[2]: volatility 1,     para[3]: volatility 2
							//               para[4]: correlation
			  )
			  // stochastic variable: R(t,T) = integral of r(u) between t and T
			  // return V(t,T) = variance of R(t,T) conditional of F_t under HW 2 Factor Model
{
	double ea, eb, eab;
	double sum1, sum2, sum3;

	ea 		= exp(-para[0]*ddt);
	eb 		= exp(-para[1]*ddt);
	eab 	= exp(-(para[0]+para[1])*ddt);

	sum1	= ddt +2.0/para[0]*ea - 0.5/para[0]*ea*ea -1.5/para[0];
	sum2	= ddt +2.0/para[1]*eb - 0.5/para[1]*eb*eb -1.5/para[1];
	sum3	= ddt +(ea-1.0)/para[0] + (eb-1.0)/para[1] -(eab-1.0)/(para[0]+para[1]);

	return	  para[2]*para[2]/para[0]/para[0]*sum1
		+ para[3]*para[3]/para[1]/para[1]*sum2
		+ 2.0*para[4]*para[2]*para[3]/para[0]/para[1]*sum3;
}

double k_PtT2(
			double t,        // start time t
			double T,        // maturity
			double dt,       // time-step for finite-difference
			double x,        // status variable x for HW 2 factor model
			double y,        // status variable y for HW 2 factor model
			                 // short rate r(t) = x(t) + y(t) + phi2(t)
			double *para,    // HW parameters para[0]: mean reversion 1, para[1]: mean reversion 2
							 //               para[2]: volatility 1,     para[3]: volatility 2
							 //               para[4]: correlation
			k_cZeroCurve *ZC  // Zero Cureve Structure
			)
// zero coupon bond price at time t maturing at T with unit face value 1 
//  under HW 2 factor model
{
	double A;


	A = 0.5* (k_Vddt2(T-t,para) - k_Vddt2(T,para) + k_Vddt2(t,para)) - (1.0-exp(-para[0]*(T-t)))/para[0]*x
		- (1.0-exp(-para[1]*(T-t)))/para[1]*y;


	return (*ZC).DF(T)/(*ZC).DF(t)*exp(A);

}

void k_FDM_backward_induction_HW1Fpp_OS_const_coeff(
	double **V_init,      // output: security values at t_init
	double **V_term,      // input: security values at t_term
	double  dt,           // time step
	                      // x: status variable for short rate(discounting zero curve)
	double  dx,           // grid step size of x
	int		Nx,           // Number of grids of x
	double  xmin,         // minimum of x
                          // y: status variable for the other short rate 
	double  dy,		      // grid step size of y
	int		Ny,           // Number of grids of y
	double  ymin,         // minimum of y
	double  *kpara,       // HW 1 factor parameters for short rate(discounting zero curve)
	                      // kpara[0]=mean reversion, kpara[1]=volatility
	double  *epara,       // HW 1 factor parameters for the other zero curve
	                      // epara[0]=mean reversion, epara[1]=volatility
	double  rho,          // correlation between two short rates
	double  *r            // short rates at t_init w.r.t x
	                      // r[i] = x[i] + phi1(t_init) 
	)
	// 2-Dim FDM backward induction for HW1(x)+HW1(y) PDE
	//           based on HangSeob's 2-Dim FDM backward induction procedure
	//
	// x: status variable for discounting curve
    // y: the other status variable
	// Quanto case에는 epara를 quanto adjustment에 의해 적절히 변형하여 투입
	//
	//	one time step backward induction
	//	find V_init with given V_term with one time step from t_term to t_init
	//
	//	Underlying PDE:
	//		V_t + c1*V_xx + c2*V_xy + c3*V_yy + c4*V_x + c5*V_y + c6*V = 0
	//		
	//	Gamma zero condition as boundary condition
	
{
	int		i, j, I, J;
	double 	*A, *B, *C, *q;
	double	**V_btwn;		//	V between t_init and t_term ; V at t = ( t_init + t_term)/2

	double c1,c2,c3,*c4,*c5;

	V_btwn	= k_misc_new_double_array_HS( Nx+1, Ny+1);

	I = Nx - 1;
	J = Ny - 1;

	// first sweep  //////////////////////////////////////////////////////////////////////////////////////////////

	A = new double [I+1];
	B = new double [I+1];
	C = new double [I+1];
	q = new double [I+1];

	c4 = new double [I+1];
	c5 = new double [J+1];

	c1 = 0.5*kpara[1]*kpara[1];
	c2 = rho*kpara[1]*epara[1];
	c3 = 0.5*epara[1]*epara[1];


	for (i=0; i <= I; i++){
		c4[i] = -kpara[0]*(xmin+i*dx);

	}

	for (j=0; j <= J; j++)
		c5[j] = -epara[0]*(ymin+j*dy);



	for ( j = 1; j <= J; j++ ) {

		for ( i = 1; i <= I; i++ ) {

			A[i] = c1/(dx*dx) - c4[i]/(2.*dx);
			B[i] = -1./dt - 2.*c1/(dx*dx) - r[i];
			C[i] = c1/(dx*dx) + c4[i]/(2.*dx);
			q[i] = 	-0.5*c2*(V_term[i+1][j+1] - V_term[i+1][j-1] - V_term[i-1][j+1] + V_term[i-1][j-1])
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

			A[j] = c3/(dy*dy) - c5[j]/(2.*dy);
			B[j] = -1./dt - 2.*c3/(dy*dy);
			C[j] = c3/(dy*dy) + c5[j]/(2.*dy);
			q[j] = 	-0.5*c2*(V_btwn[i+1][j+1] - V_btwn[i+1][j-1] - V_btwn[i-1][j+1] + V_btwn[i-1][j-1])
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
	if (c4) delete[] c4;
	if (c5) delete[] c5;


	k_misc_delete_array_HS( V_btwn, Nx+1 );
}
