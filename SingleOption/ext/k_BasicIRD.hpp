#ifndef KDB_BASICIRD_H
#define KDB_BASICIRD_H

#include "kdb_date.hpp"
#include "kdb_swap_exact_schedule.hpp"

/*******************************************************************************************/
/* Program Name : kdb_BasicIRD.hpp
/* Version : 0.0.0.0
/* Author : Great Quant Jaemin Ahn
/* Date : Sep. 1, 2010
/* Modified by : 
/* Modified at :  
/* Copyright : KDB Quantitative Analysis Team 
/* Description : A collection of basic interest rate derivatives and related functions
/********************************************************************************************/




double k_Black_d1( double F_0, double K, double T, double sigma );
    // parameter for Black-Scholes Formula  
	//	F_0	:	forward price (maturity T, seen at t=0)
	//	K	:	strike
	//	T	:	option maturity 
	//	sigma	:	volatility of underlying


double k_BlackCall( double disBond_0_Tstar, double F_0, double K, double T, double sigma);
    // Black-Scholes Call Formula 
	//	Tstar	:	payoff delivery time
	//	disBond_0_Tstar	:	zero-coupon bond price maturing at Tstar
	//	F_0	:	forward price (maturity T, seen at t=0)
	//	K	:	strike
	//	T	:	option maturity
	//	sigma	:	volatility of underlying 


double k_BlackPut( double disBond_0_Tstar, double F_0, double K, double T, double sigma);
    // Black-Scholes Put Formula 
	//	Tstar	:	payoff delivery time
	//	disBond_0_Tstar	:	zero-coupon bond price maturing at Tstar
	//	F_0	:	forward price (maturity T, seen at t=0)
	//	K	:	strike
	//	T	:	option maturity
	//	sigma	:	volatility of underlying

double k_caplet( double disBond_mat, double F, double K, double c_s, double c_int, double L, double sigma);
	//	disBond_mat	:	disFactor from today to payment date(calculation end date)
	//                     or Zero Coupon Bond price maturing at payment date (face value 1)
	//  F           :   today's forward rate 
	//	K			:	strike
	//	c_s			:	calculation start date( caplet maturity)
	//	c_int	    :	calculation period size = calculation end - calculation start 
	//                                          = payment date - caplet maturity
	//                                              considering day-count convention
	//	L			:	principal
	//	sigma		:	volatility


double k_floorlet( double disBond_mat, double F, double K, double c_s, double c_int, double L, double sigma);
//	disBond_mat	:	disFactor from today to payment date(calculation end date)
//                     or Zero Coupon Bond price maturing at payment date (face value 1)
//  F           :   today's forward rate 
//	K			:	strike
//	c_s			:	calculation start date( floorlet maturity)
//	c_int	    :	calculation period size = calculation end - calculation start 
//                                          = payment date - floorlet maturity
//                                              considering day-count convention
//	L			:	principal
//	sigma		:	volatility

double k_cap(
			 k_cZeroCurve *zCurve,
			 double strike, 
			 double L,
			 double sigma,
			 date today,
			 date capmat,
			 PayType ptype,
			 int BizDayConv,
			 int DayCountConv,
			 calendar pay_cal,
			 calendar fix_cal
			 );
	//	return price of cap using flat vol
	//	the caplet whose payoff fixed at today is not included

	//	*zCurve : zero curve structure, 	
	//	strike : cap strike 
	//	L : principal
	//  sigma : cap volatility
	//	today :	today (trade date)
	//  capmat  : cap maturity 
	//  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually
	//  BizDayConv : 0=modified following,  1=following,  2=preceding
	//	DayCountConv :  0=A/365, 1=30/360, 2=A/360, 3=A/A
	//  pay_cal : payment calendar --> SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, 
	//                                 SEOUL_NEW_YORK, TARGET_NEW_YORK, TOKYO_NEW_YORK
	//  fix_cal : fixing calendar -->  SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, 
	//                                 SEOUL_NEW_YORK, TARGET_NEW_YORK, TOKYO_NEW_YORK
	// 현재 거래시작 시점에서의 가격만 산출

double k_cap( k_cZeroCurve *zCurve, double strike, int nb_fixing, double L, double sigma, PayType ptype );
	// simplified cap pricer 
	//	return price of cap using flat vol
	//	the caplet whose payoff fixed at today is not included

	//	nb_fixing : number of fixings = number of caplets, 	
	//	L : principal
	//	sigma: volatility
	//  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually

double k_floor(
			 k_cZeroCurve *zCurve,
			 double strike, 
			 double L,
			 double sigma,
			 date today,
			 date floormat,
			 PayType ptype,
			 int BizDayConv,
			 int DayCountConv,
			 calendar pay_cal,
			 calendar fix_cal
			 );
//	return price of cap using flat vol
//	the floorlet whose payoff fixed at today is not included

//	*zCurve : zero curve structure, 	
//	strike : floor strike 
//	L : principal
//  sigma : floor volatility
//	today :	today (trade date)
//  floormat  : floor maturity 
//  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually
//  BizDayConv : 0=modified following,  1=following,  2=preceding
//	DayCountConv :  0=A/365, 1=30/360, 2=A/360, 3=A/A
//  pay_cal : payment calendar --> SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, 
//                                 SEOUL_NEW_YORK, TARGET_NEW_YORK, TOKYO_NEW_YORK
//  fix_cal : fixing calendar -->  SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, 
//                                 SEOUL_NEW_YORK, TARGET_NEW_YORK, TOKYO_NEW_YORK
// 현재 거래시작 시점에서의 가격만 산출


double k_floor( k_cZeroCurve *zCurve, double strike, int nb_fixing, double L, double sigma, PayType ptype );
	// simplified floor pricer 
	//	return price of floor using flat vol
	//	the floorlet whose payoff fixed at today is not included

	//	nb_fixing : number of fixings = number of floorlets, 	
	//	L : principal
	//	sigma: volatility
	//  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually

double k_ForwardSwapRate(
						k_cZeroCurve *zCurve,
						date today,
						date omat,
						date smat,
						PayType ptype,
						FRateType frtype,
						int BizDayConv,
						int DayCountConv,
						calendar pay_cal,
						calendar fix_cal
						);
	//  No specific short rate model is not involved.
	//  just generate forward swap rate and spot swap rate from zero curve
	//  !!caution:  "kdb_date.hpp" && "kdb_swap_exact_schedule.hpp" should be included
	//  *zCurve : zero curve structure
	//	today :	today
	//  omat  : swap trade date
    //  smat  : underlying swap maturity 
	//  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually
	//  frtype : Floating Rate Type  M3=3 Month, M6=6 Month
	//  BizDayConv : 0=modified following,  1=following,  2=preceding
	//	DayCountConv :  0=A/365, 1=30/360, 2=A/360, 3=A/A
	//  pay_cal : payment calendar --> SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, SEOUL_NEW_YORK,
	//                                 TARGET_NEW_YORK, TOKYO_NEW_YORK
	//  fix_cal : fixing calendar --> SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, SEOUL_NEW_YORK,
	//                                TARGET_NEW_YORK, TOKYO_NEW_YORK


double k_ForwardSwapRate( int n, double *Tau, double *ZB);
	//	generate forward swap rate and spot swap rate from array of zero coupon bond price
	//	formular derived by hangseob	:	need to check!!!	
	//	n:	number of payment
	//	Tau[0] ~ Tau[n-1]: payment interval size considering day count convention
	//                     Tau[0] = first payment day - swap starting day
	//                     Tau[1] = second payment day - first payment day
	//                 ... Tau[n-1] = n-th payment day - (n-1)-th payment day 
	//	ZB[k]: zero coupon bond price maturing at the (k + 1)-th payment day


double k_ForwardSwapRate( int n, double tau, double *ZB);
	//	generate forward swap rate and spot swap rate from array of zero coupon bond price
	//	formular derived by hangseob	:	need to check!!!
	//  
	//	n:	number of payment
	//	tau: fixed payment interval, namely 0.25 = 3M, 0.5 = 6M, ....
	//	ZB[k]: zero coupon bond price maturing at the (k + 1)-th payment day


double k_ForwardSwapRate(k_cZeroCurve *zCurve, int n, double Ts, PayType ptype);
	// simplified version;
	//	generate forward swap rate and spot swap rate from zero curve
	//	formular derived by hangseob	:	need to check!!!	modified by Jaemin Ahn
	//	n:	number of payment
	//	Ts : swap starting time(i.e. first fixing point)
	//		*caution!!: when it is fwd swap, Ts is not zero!!!
	//  ptype : Payment Type(fixed coupon) M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually


double k_Swaption( 
				  int Flag,
				  double strike,
				  double sigma,
				  double L,
				  k_cZeroCurve *zCurve,
				  date today,
				  date omat,
				  date smat,
				  PayType ptype,
				  FRateType frtype,
				  int BizDayConv,
				  int DayCountConv,
				  calendar pay_cal,
				  calendar fix_cal
				  )  ;
	//	if option holder pays fixed, Flag = 0;
	//	if option holder receives fixed, Flag = 1;
	//	strike :	strike of swaption
	//  sigma : volatility 
	//  L : Notional amount
	//  *zCurve : zero curve structure
	//	today :	today
	//  omat  : swap trade date
	//  smat  : underlying swap maturity 
	//  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually
	//  frtype : Floating Rate Type  M3=3 Month, M6=6 Month
	//  BizDayConv : 0=modified following,  1=following,  2=preceding
	//	DayCountConv :  0=A/365, 1=30/360, 2=A/360, 3=A/A
	//  pay_cal : payment calendar --> SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, SEOUL_NEW_YORK,
	//                                 TARGET_NEW_YORK, TOKYO_NEW_YORK
	//  fix_cal : fixing calendar --> SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, SEOUL_NEW_YORK,
	//                                TARGET_NEW_YORK, TOKYO_NEW_YORK


double k_Swaption(
				  int Flag,
				  k_cZeroCurve *zCurve,
				  int n,
				  double Ts,
				  PayType	ptype,
				  double strike,
				  double sigma,
				  double L)  ;
	//   simplified swaption pricer
	//	if option holder pays fixed, Flag = 0;
	//	if option holder receives fixed, Flag = 1;
	//	n:	number of payment
	//	Ts :	swap starting time(i.e option maturity time)
	//   ptype: fixed rae payement type M=1M, Q=3M, S=6M, A=12M
	//   strike: swaption strike
	//   sigma : volatility
	//   L : principal


double k_caplet_stripping(
						  k_cZeroCurve *zCurve,
						  double strike,
						  date today,
						  double *mk_cap_vol,
						  date *mk_cap_mat,
						  int num_cap,
						  PayType ptype,
						  int BizDayConv,
						  int DayCountConv,
						  calendar pay_cal,
						  calendar fix_cal,
						  int num_caplet,
						  date *caplet_mat,
						  double *caplet_vol
						  );

  //  *zCurve : zero curve structure, 	
  //  strike : cap strike 
  //  today :	today (trade date)
  //  *mk_cap_vol : market cap volatilities w.r.t. maturities 
  //  *mk_cap_mat : cap maturities
  //  num_cap : number of caps (number of maturities) 
  //  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually
  //  BizDayConv : 0=modified following,  1=following,  2=preceding
  //	DayCountConv :  0=A/365, 1=30/360, 2=A/360, 3=A/A
  //  pay_cal : payment calendar --> SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, 
  //                                 SEOUL_NEW_YORK, TARGET_NEW_YORK, TOKYO_NEW_YORK
  //  fix_cal : fixing calendar -->  SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, 
  //                                 SEOUL_NEW_YORK, TARGET_NEW_YORK, TOKYO_NEW_YORK
  //  num_caplet = number of caplets 
  //  Output: *caplet_mat = caplet maturities 
  //          *caplet_vol = stripped caplet volatilities 


double k_caplet_stripping(
						  k_cZeroCurve *zCurve,
						  double strike,
						  double *mk_cap_vol,
						  double *mk_cap_mat,
						  int num_cap,
						  PayType ptype,
						  int num_caplet,
						  double *caplet_mat,
						  double *caplet_vol
						  );
  // simplified caplet_stripping 
  //	the caplet whose payoff fixed at today is not included
  //	*zCurve : zero curve structure, 	
  //	strike : cap strike 
  //  *mk_cap_vol : market cap volatilities w.r.t. maturities 
  //  *mk_cap_mat : cap maturities
  //  num_cap : number of caps (number of maturities) 
  //  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually
  //  num_caplet : number of caplets 
  //   Output:*caplet_mat = caplet maturities
  //          *caplet_vol = stripped caplet volatilities 



#endif