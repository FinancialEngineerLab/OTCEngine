/*******************************************************************************************/
/* Program Name : kdb_HW_functions.hpp
/* Version : 1.0.0.0
/* Author : Great Quant Jaemin Ahn
/* Date : Jan. 12, 2011
/* Modified by : 
/* Modified at :  
/* Copyright : KDB Quantitative Analysis Team 
/* Description : A collection of functions related to HW model 
/********************************************************************************************/

#ifndef KDB_HW_FUNCTIONS_H
#define KDB_HW_FUNCTIONS_H

#include <kdb_intRateTermStructure.hpp>
#include <kdb_miscellaneous.hpp>
#include <kdb_LinearAlgebra.hpp>

double k_phi1( 
			double T,           // time
			double dt,          // time-step for finite-difference
			double *para,       // HW parameters para[0]: mean reversion, para[1]: volatility
			k_cZeroCurve *ZC    // Zero Curve Structure
			);
			// HW 1 Factor short rate mean level function r(T)t = x(T) + phi1(T)

double k_PtT1(
			double t,         // start time t
			double T,         // maturity
			double dt,        // time-step for finite-difference
			double r,         // short rate at t
			double *para,     // HW paraemters para[0]: mean reversion, para[1]: volatility
			k_cZeroCurve *ZC    // Zero Cureve Structure
			);
			// zero coupon bond price at time t maturing at T with unit face value 1 
			//  under HW 1 factor model



double	k_Vddt1( 
			  double ddt,   // T-t
			  double *para  // HW parameters para[0]: mean reversion, para[1]: volatility 
			  );
			  // stochastic variable: R(t,T) = integral of r(u) between t and T
              // return V(t,T) = variance of R(t,T) conditional of F_t under HW 1 Factor model


double	k_phi2( 
			 double T,           // time
			 double dt,          // time-step for finite-difference
			 double *para,       // HW parameters para[0]: mean reversion 1, para[1]: mean reversion 2
			 //               para[2]: volatility 1,     para[3]: volatility 2
			 //               para[4]: correlation
			 k_cZeroCurve *ZC    // Zero Curve Structure
			 );
			 // HW 2 Factor short rate mean level function r(T) = x(T) + y(T) + phi2(T)


double	k_Vddt2(
			  double ddt,   // T-t
			  double *para  // HW paraemters para[0]: mean reversion 1, para[1]: mean reversion 2
			  //               para[2]: volatility 1,     para[3]: volatility 2
			  //               para[4]: correlation
			  );
			 // stochastic variable: R(t,T) = integral of r(u) between t and T
			 // return V(t,T) = variance of R(t,T) conditional of F_t under HW 2 Factor Model

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
			);
			// zero coupon bond price at time t maturing at T with unit face value 1 
			//  under HW 2 factor model


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
	double  *kpara,       // HW 1 factor parameters for the discounting zero curve
	// kpara[0]=mean reversion, kpara[1]=volatility
	double  *epara,       // HW 1 factor parameters for the other zero curve
	// epara[0]=mean reversion, epara[1]=volatility
	double  rho,          // correlation between two short rates
	double  *r            // short rates at t_init w.r.t x
	// r[i] = x[i] + phi1(t_init) 
	);
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
	
#endif