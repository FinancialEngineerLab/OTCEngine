#ifndef KDB_OPTION_FORMULA_HPP
#define KDB_OPTION_FORMULA_HPP

/*************************************************************************
*	Program Name : kdb_OptionFormula.hpp
* Version : 0.0.0.2
* Author : Hwang-hyun Kwon
* Date : Jan. 15, 2009
* Last Modified by : Great Quant Jaemin Ahn
* Last Modified at : Feb. 22, 2011 
* Copyright : KDB Quantitative Analysis Team 
* Description : A collection of closed-form option formulae
**************************************************************************/

double k_fmla_BSCall( double dSpot,
					 double dStrike,
					 double dRate,
					 double dDividend,
					 double dVol,
					 double dExpiry);

double k_fmla_BSPut( double dSpot,
					double StRateike,
					double dRate,
					double dDividend,
					double dVol,
					double dExpiry);


double k_fmla_BSDigitalCall(double dSpot,
							double dStrike,
							double dRate,
							double dDividend,
							double dVol,
							double dExpiry);

double k_fmla_BSDigitalPut(double dSpot,
						   double dStrike,
						   double dRate,
						   double dDividend,
						   double dVol,
						   double dExpiry);

double k_fmla_BSCallDelta( double dSpot,
						  double dStrike,
						  double dRate,
						  double dDividend,
						  double dVol,
						  double dExpiry);

double k_fmla_BSPutDelta( double dSpot,
						 double dStrike,
						 double dRate,
						 double dDividend,
						 double dVol,
						 double dExpiry);						 

double k_fmla_BSGamma( double dSpot,
					  double dStrike,
					  double dRate,
					  double dDividend,
					  double dVol,
					  double dExpiry);


double k_fmla_BSVega( double dSpot,
					 double dStrike,
					 double dRate,
					 double dDividend,
					 double dVol,
					 double dExpiry);

double k_fmla_BSCallTheta( double dSpot,
						  double dStrike,
						  double dRate,
						  double dDividend,
						  double dVol,
						  double dExpiry);

double k_fmla_BSPutTheta(double dSpot, double dStrike, double dRate, double dDividend, double dVol, double dExpiry) ; 

double k_fmla_BSCallRho( double dSpot,
						double dStrike,
						double dRate,
						double dDividend,
						double dVol,
						double dExpiry);

double k_fmla_BSPutRho(double dSpot, double dStrike, double dRate, double dDividend, double dVol, double dExpiry) ; 	

double k_fmla_BlackCall(double dForward, double dStrike, double dRate, double dVol, double dExpiry); 

double k_fmla_BlackPut(double dForward, double dStrike, double dRate, double dVol, double dExpiry); 

double k_fmla_CallOnCall(double dSpot, double dStrike, double dStrikeOp, double dRate, double dDividend, 
						 double dVol, double dExpiry, double dExpiryOp) ; 

double k_fmla_GAvg(int iCallOrPut, double dSpot, double dStrike, double dRate, double dDividend, double dVol, double dExpiry) ; 

double k_fmla_AAvgLevy(int iCallOrPut, double dSpot, double dSpotAvg, double dStrike, double dRate, double dCostCarry, 
					   double dVol, double dExpiry, double dRemainingExpiry) ; 

double k_fmla_SpreadAprox(int iCallOrPut, double dSpot1, double dSpot2, double dAmount1, double dAmount2, double dCorr, double dStrike, 
						  double dRate, double dCostCarry1, double dCostCarry2, double dVol1, double dVol2, double dExpiry) ; 

void set_intermediate_variables(double S, double X, double H, double r, double div, double sigma, double T, double K,		
								double eta, double phi,																					
								double *Acap, double *Bcap, double *Ccap, double *Dcap, double *Ecap, double *Fcap);

double k_fmla_Barrier( long call_or_put, long up_or_down, long in_or_out,  
					  double S, double X, double H, double r, double div, double sigma, double T, double K);




double k_KI_CashOrNothing(double S, double H, double r, double div, double sigma, double T, double K, int UD);
	// Knock-In Cash (at expiration) or nothing option 
	// S: spot
	// H: barrier
	// r: risk free interest rate
	// div: continuous dividend rate
	// sigma: volatility
	// T: time to expiration
	// K: rebate 	
	// UD = 1: Up,  UD = 0: Down


double k_KO_Rebate_at_Expiry(double S, double K, double H, double r, double div, double sigma, double T, double R, int UD, int CP);
	// KO option with rebate at expiry 


// Double Barrier Global KO option
double k_DoubleBarrierKO(
						 double S,            // Spot
						 double K,            // Strike
						 double U,            // Upper Barrier
						 double L,            // Lower Barrier
						 double T,            // time to maturity
						 double r,            // risk-free rate
						 double q,            // dividend rate
						 double vol,          // volatility
						 int CP               // CP = 1: Call    CP=0: Put
						 ) ;
						 
// Double Barrier Global KI option
double k_DoubleBarrierKI(
						 double S,            // Spot
						 double K,            // Strike
						 double U,            // Upper Barrier
						 double L,            // Lower Barrier
						 double T,            // time to maturity
						 double r,            // risk-free rate
						 double q,            // dividend rate
						 double vol,          // volatility
						 int CP               // CP = 1: Call    CP=0: Put
						 );

// One Touch Double Barrier KI (rebate R at maturity)
double k_OneTouchDoubleBarrierKI(
								 double S,            // Spot
								 double U,            // Upper Barrier
								 double L,            // Lower Barrier
								 double T,            // time to maturity
								 double r,            // risk-free rate
								 double q,            // dividend rate
								 double vol,          // volatility
								 double R             // Rebate
								 );

// One Touch Double Barrier KO (rebate R at maturity)
double k_OneTouchDoubleBarrierKO(
								 double S,            // Spot
								 double U,            // Upper Barrier
								 double L,            // Lower Barrier
								 double T,            // time to maturity
								 double r,            // risk-free rate
								 double q,            // dividend rate
								 double vol,          // volatility
								 double R             // Rebate
								 );

double integrand_for_spread_option(double z, double K, double x1, double x2, double x11, double x22, double x12,
								   double b1, double b2, double nu, int CP);
								   // integrand for spread option approximation

double adaptive_quad_for_spread_option(double a, double b, double fa, double fb, double TOL,
									   double K, double x1, double x2, double x11, double x22, double x12,
									   double b1, double b2, double nu, int *level,int CP);
									   // adaptive quadrature for spread option approximation

double k_Spread_Option(
					   double Disfactor, // discount factor from spot day to maturity(delivery date)
					   double F1, // today's forward price for S1 
					   double F2,  // today's forward price for S2 						  
					   double K,    // Strike
					   double v1,   // volatility for S1
					   double v2,   // volatility for S2
					   double rho,  // correlation
					   double r,    // risk-free rate
					   double T,     // maturity
					   int CP        // CP = 1: Call,  CP = 0: Put
						  );
    // Spread Option Approximation based on MUREX spread option formula
double call_iv(double dSpot,
					 double dStrike,
					 double dRate,
					 double dDividend,
					 double dVol,
					 double dExpiry,
					 double dCallPrice);
double put_iv(double dSpot,
					 double dStrike,
					 double dRate,
					 double dDividend,
					 double dVol,
					 double dExpiry,
					 double dPutPrice);
#endif

