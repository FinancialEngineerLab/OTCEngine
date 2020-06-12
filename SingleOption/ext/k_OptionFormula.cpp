
#include <cmath>
#include <algorithm>
#include "k_OptionFormula.hpp"
#include "k_probability.hpp"
#include "k_BasicMath.hpp"
#include "assert.h"



using namespace std; 

double k_fmla_BSCall( double dSpot,
					 double dStrike,
					 double dRate,
					 double dDividend,
					 double dVol,
					 double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);

	if (standardDeviation == 0.0){
		return(exp(-dRate*dExpiry)*max(dSpot-dStrike,0.0) );
	}
	double moneyness = log(dSpot/dStrike);
	double d1 =( moneyness +  (dRate-dDividend)*dExpiry+0.5* standardDeviation*standardDeviation)/standardDeviation;
	double d2 = d1 - standardDeviation;
	return dSpot*exp(-dDividend*dExpiry) * k_pr_CN(d1) - dStrike*exp(-dRate*dExpiry)*k_pr_CN(d2);
}
// Black Scholes Call Price 
// the meanings of the variables should be obvious from their names. 
// (risk free) rate, dividend, vol are all in real numbers, not in percentage. 
// dexpiry stands for "time to expiry" in years, e.g. 1.25 for 1 yr and 3 month.  
// and the same for the following functions. 

double k_fmla_BSPut( double dSpot,
					double dStrike,
					double dRate,
					double dDividend,
					double dVol,
					double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	
	if (standardDeviation == 0.0){
		return(exp(-dRate*dExpiry)*max(dStrike-dSpot,0.0) );
	}

	//return(exp(-dRate*dExpiry)*max(dSpot-dStrike,0.0) );
	double moneyness = log(dSpot/dStrike);
	double d1 =( moneyness +  (dRate-dDividend)*dExpiry+0.5* standardDeviation*standardDeviation)/standardDeviation;
	double d2 = d1 - standardDeviation;
	return dStrike*exp(-dRate*dExpiry)*(1.0-k_pr_CN(d2)) - dSpot*exp(-dDividend*dExpiry) * (1-k_pr_CN(d1));
}


// Black Scholes formulae for digital options
double k_fmla_BSDigitalCall( double dSpot,
							double dStrike,
							double dRate,
							double dDividend,
							double dVol,
							double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);

	if (standardDeviation == 0.0){
		if (dSpot >= dStrike)
			return(exp(-dRate*dExpiry));
		else
			return(0.0);
	}

	double moneyness = log(dSpot/dStrike);
	double d2 =( moneyness +  (dRate-dDividend)*dExpiry-0.5* standardDeviation*standardDeviation)/standardDeviation;
	return exp(-dRate*dExpiry)*k_pr_CN(d2);
}
// also called "cash or nothing" call option
// paying one unit of cash if the spot is above the strike at maturity

double k_fmla_BSDigitalPut( double dSpot,
						   double dStrike,
						   double dRate,
						   double dDividend,
						   double dVol,
						   double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);

	if (standardDeviation == 0.0){
		if (dSpot <= dStrike)
			return(exp(-dRate*dExpiry));
		else
			return(0.0);
	}

	double moneyness = log(dSpot/dStrike);
	double d2 =( moneyness +  (dRate-dDividend)*dExpiry-0.5* standardDeviation*standardDeviation)/standardDeviation;
	return exp(-dRate*dExpiry)*(1.0-k_pr_CN(d2));
}


// Black Scholes Greeks
double k_fmla_BSCallDelta( double dSpot,
						  double dStrike,
						  double dRate,
						  double dDividend,
						  double dVol,
						  double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dSpot/dStrike);
	double d1 =( moneyness +  (dRate-dDividend)*dExpiry+0.5* standardDeviation*standardDeviation)/standardDeviation;
	double result = exp(-dDividend*dExpiry)*k_pr_CN(d1) ; 
	return result; 
}								 

double k_fmla_BSPutDelta( double dSpot,
						 double dStrike,
						 double dRate,
						 double dDividend,
						 double dVol,
						 double dExpiry)					 
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dSpot/dStrike);
	double d1 =( moneyness +  (dRate-dDividend)*dExpiry+0.5* standardDeviation*standardDeviation)/standardDeviation;
	double result = exp(-dDividend*dExpiry)*(k_pr_CN(d1) -1.0); 
	return result; 
}
double k_fmla_BSGamma( double dSpot,
					  double dStrike,
					  double dRate,
					  double dDividend,
					  double dVol,
					  double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dSpot/dStrike);
	double d1 =( moneyness +  (dRate-dDividend)*dExpiry+0.5* standardDeviation*standardDeviation)/standardDeviation;
	double result = exp(-dDividend*dExpiry)*k_pr_Npdf(0.0,1.0,d1)/(dSpot*standardDeviation) ; 
	return result; 
}	

double k_fmla_BSVega( double dSpot,
					 double dStrike,
					 double dRate,
					 double dDividend,
					 double dVol,
					 double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dSpot/dStrike);
	double d1 =( moneyness +  (dRate-dDividend)*dExpiry+0.5* standardDeviation*standardDeviation)/standardDeviation; 
	double NormalDensity = ONE_OVER_SQRT_TWO_PI *exp(-d1*d1/2);
	return dSpot*exp(-dDividend*dExpiry) * sqrt(dExpiry)*NormalDensity ;
}


double k_fmla_BSCallTheta(double dSpot, double dStrike, double dRate, double dDividend, double dVol, double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dSpot/dStrike);
	double d1 =( moneyness +  (dRate-dDividend)*dExpiry+0.5* standardDeviation*standardDeviation)/standardDeviation;
	double d2 = d1 - standardDeviation;
	double result = -dSpot*exp(-dDividend*dExpiry)*(0.5*k_pr_Npdf(0.0, 1.0, d1)*dVol/sqrt(dExpiry) - dDividend*k_pr_CN(d1) )
		-dRate*dStrike*exp(-dRate*dExpiry)*k_pr_CN(d2) ; 
	return result; 
}


double k_fmla_BSPutTheta(double dSpot, double dStrike, double dRate, double dDividend, double dVol, double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dSpot/dStrike);
	double d1 =( moneyness +  (dRate-dDividend)*dExpiry+0.5* standardDeviation*standardDeviation)/standardDeviation;
	double d2 = d1 - standardDeviation;
	double result = -dSpot*exp(-dDividend*dExpiry)*(0.5*k_pr_Npdf(0.,1.,d1)*dVol/sqrt(dExpiry) + dDividend*k_pr_CN(-d1) )
		+dRate*dStrike*exp(-dRate*dExpiry)*k_pr_CN(-d2) ; 
	return result; 
}	


double k_fmla_BSCallRho(double dSpot, double dStrike, double dRate, double dDividend, double dVol, double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dSpot/dStrike);
	double d2 =( moneyness +  (dRate-dDividend)*dExpiry-0.5* standardDeviation*standardDeviation)/standardDeviation;
	double result = dExpiry*dStrike*exp(-dRate*dExpiry)*k_pr_CN(d2) ; 
	return result; 
}


double k_fmla_BSPutRho(double dSpot, double dStrike, double dRate, double dDividend, double dVol, double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dSpot/dStrike);
	double d2 =( moneyness +  (dRate-dDividend)*dExpiry-0.5* standardDeviation*standardDeviation)/standardDeviation;
	double result = -dExpiry*dStrike*exp(-dRate*dExpiry)*k_pr_CN(-d2) ; 
	return result; 
}	


// Black formulae for options on futures or forwards 
double k_fmla_BlackCall(double dForward, double dStrike, double dRate, double dVol, double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dForward/dStrike);
	double d1 =( moneyness +0.5* standardDeviation*standardDeviation)/standardDeviation;
	double d2 = d1 - standardDeviation;
	double result = exp(-dRate*dExpiry) *(dForward*k_pr_CN(d1) - dStrike*k_pr_CN(d2));
	return result; 
}


double k_fmla_BlackPut(double dForward, double dStrike, double dRate, double dVol, double dExpiry)
{
	double standardDeviation = dVol*sqrt(dExpiry);
	double moneyness = log(dForward/dStrike);
	double d1 =( moneyness +0.5* standardDeviation*standardDeviation)/standardDeviation;
	double d2 = d1 - standardDeviation;
	double result = exp(-dRate*dExpiry) *(dStrike*k_pr_CN(-d2) - dForward*k_pr_CN(-d1));
	return result; 
}


// "Approximate" Formulae for the Call on Call Option using Black Scholes ones
// based on the work by Bensoussan, Crouhy, and Galai('95). see p.136 of Haug's Option formula book.  
// "dExpiryOp" and "dStrikeOp" stand for the time to expiry and the strike of the option on option, respectively.  
double k_fmla_CallOnCall(double dSpot, double dStrike, double dStrikeOp, double dRate, double dDividend, double dVol, double dExpiry, double dExpiryOp)
{
	double dCallValue = k_fmla_BSCall(dSpot, dStrike, dRate,  dDividend,  dVol,  dExpiry) ; 
	double dDelta = k_fmla_BSCallDelta( dSpot,  dStrike,  dRate,  dDividend,  dVol,  dExpiry) ; 
	double dAbsDelta = fabs(dDelta) ; 
	double dOptionVol = dVol*dAbsDelta*dSpot/dCallValue ; 
	double dStdDev = dOptionVol*sqrt(dExpiryOp) ; 
	double d1 = (log(dCallValue/dStrikeOp) + (dRate-dDividend)*dExpiryOp +0.5*dStdDev*dStdDev)/dStdDev ;
	double d2 = d1 - dStdDev ; 
	double result = dCallValue*k_pr_CN(d1) - dStrikeOp*exp(-dRate*dExpiryOp)*k_pr_CN(d2) ; 
	return result; 
}	

// Asian Options 
// Geometric Continuous Average-Rate Call and Put Options
// iCallOrPut = 1 for call and -1 for put 
double k_fmla_GAvg(int iCallOrPut, double dSpot, double dStrike, double dRate, double dDividend, double dVol, double dExpiry) 
{
	int isign = iCallOrPut ; 
	double dAdjVol = dVol/sqrt(3.0) ; 
	double dCostCarry = 0.5*(dRate-dDividend - dVol*dVol/6.0) ; 
	double dStdDev = dAdjVol*sqrt(dExpiry) ; 
	double d1 = (log(dSpot/dStrike) + dCostCarry*dExpiry + 0.5*dStdDev*dStdDev)/dStdDev ;
	double d2 = d1 -dStdDev ; 
	double result = isign*(dSpot*exp((dCostCarry-dRate)*dExpiry)*k_pr_CN(isign*d1) 
		- dStrike*exp(-dRate*dExpiry)*k_pr_CN(isign*d2) ) ; 
	return result ; 
}

// Levy Approximation for Arithmetic Asian Option
// iCallOrPut = 1 for call and -1 for put 
double k_fmla_AAvgLevy(int iCallOrPut, double dSpot, double dSpotAvg, double dStrike, double dRate, double dCostCarry, 
					   double dVol, double dExpiry, double dRemainingExpiry)
					   // dSpotAvg = arithmetic average of the known asset price fixings
					   // dCostCarry = cost-of-carry rate, (risk-free rate - dividend) for equities and (domestic rate - foreign rate) for FX, e.g. 
					   // dRemainingExpiry = remaing time to expiry			
					   //******** WARNING : NOT APPLICABLE WHEN COST-OF-CARRY RATE IS ZERO ***********											
{
	assert( dCostCarry != 0.0 ) ; // the formula not applicable for this case 
	int isign = iCallOrPut; 
	double dModSpot = dSpot*(exp((dCostCarry-dRate)*dRemainingExpiry) - exp(-dRate*dRemainingExpiry) )/(dExpiry*dCostCarry) ; 
	double dModStrike = dStrike - (dExpiry-dRemainingExpiry)*dSpotAvg/dExpiry ; 
	double dM = 2*dSpot*dSpot/(dCostCarry+dVol*dVol)*( (exp((2*dCostCarry+dVol*dVol)*dRemainingExpiry) -1)
		/(2*dCostCarry +dVol*dVol)  - (exp(dCostCarry*dRemainingExpiry) -1)/dCostCarry )/(dExpiry*dExpiry) ; 
	double dModVariance = log(dM) -2*(dRate*dRemainingExpiry +log(dModSpot) ) ; 
	double d1 = (0.5*log(dM) -log(dModStrike))/sqrt(dModVariance) ; 
	double d2 = d1 - sqrt(dModVariance) ; 
	double result = isign*(dModSpot*k_pr_CN(isign*d1) 
		- dModStrike*exp(-dRate*dRemainingExpiry)*k_pr_CN(isign*d2) ); 
	return result ; 								
}												


// Kirk Approximation for Spread Option
// Haug, "Option Formula", p.213
double k_fmla_SpreadAprox(int iCallOrPut, double dSpot1, double dSpot2, double dAmount1, double dAmount2, double dCorr, double dStrike, 
						  double dRate, double dCostCarry1, double dCostCarry2, double dVol1, double dVol2, double dExpiry)
{
	int isign = iCallOrPut; 
	double dA = dAmount2*dSpot2*exp((dCostCarry2-dRate)*dExpiry) ; 
	double dB = dA + dStrike*exp(-dRate*dExpiry) ; 
	double dModSpot = dAmount1*dSpot1*exp((dCostCarry1-dRate)*dExpiry) / dB ; 
	double dF = dA/dB ; 
	double dModVol = sqrt(dVol1*dVol1 + pow(dVol2*dF, 2.0) -2*dCorr*dVol1*dVol2*dF) ; 
	double dStdDev = dModVol*sqrt(dExpiry) ; 
	double d1 = (log(dModSpot) + 0.5*dStdDev*dStdDev)/dStdDev ;
	double d2 = d1 - dStdDev ; 
	double result = isign*dB*(dModSpot*k_pr_CN(isign*d1) - k_pr_CN(isign*d2) ) ; 
	return result; 
}	


// Barrier Options


void set_intermediate_variables( 
								double S, double X, double H, double r, double div, double sigma, double T, double K,		//inputs
								double eta, double phi,																		//inputs			
								double *Acap, double *Bcap, double *Ccap, double *Dcap, double *Ecap, double *Fcap			//outputs
								)
{
	// S: spot
	// X: strike
	// H: barrier
	// r: risk free interest rate
	// div: continuous dividend rate
	// sigma: volatility
	// T: time to expiration
	// K: rebate 	

	double b, mu, lambda;
	double x1, x2, y1, y2, z;

	b = r - div;

	mu = (b - 0.5*sigma*sigma)/(sigma*sigma);

	lambda = sqrt(mu*mu + 2.0*r/(sigma*sigma));

	x1 = log(S/X)/(sigma*sqrt(T)) + (1.0 + mu)*sigma*sqrt(T);

	x2 = log(S/H)/(sigma*sqrt(T)) + (1.0 + mu)*sigma*sqrt(T);

	y1 =log(H*H/(S*X))/(sigma*sqrt(T)) + (1.0 + mu)*sigma*sqrt(T);

	y2 = log(H/S)/(sigma*sqrt(T)) + (1.0 + mu)*sigma*sqrt(T);

	z = log(H/S)/(sigma*sqrt(T)) + lambda*sigma*sqrt(T);

	*Acap = phi*S*exp((b-r)*T)*k_pr_CN(phi*x1)

		- phi*X*exp(-r*T)*k_pr_CN(phi*x1 - phi*sigma*sqrt(T));

	*Bcap = phi*S*exp((b-r)*T)*k_pr_CN(phi*x2)

		- phi*X*exp(-r*T)*k_pr_CN(phi*x2 - phi*sigma*sqrt(T));

	*Ccap = phi*S*exp((b-r)*T)*pow(H/S, 2.0*(mu+1.0))*k_pr_CN(eta*y1)

		- phi*X*exp(-r*T)*pow(H/S,2.0*mu)*k_pr_CN(eta*y1 -eta*sigma*sqrt(T));	

	*Dcap = phi*S*exp((b-r)*T)*pow(H/S, 2.0*(mu+1.0))*k_pr_CN(eta*y2)

		- phi*X*exp(-r*T)*pow(H/S,2.0*mu)*k_pr_CN(eta*y2 -eta*sigma*sqrt(T));

	*Ecap = K*exp(-r*T)*(k_pr_CN(eta*x2 - eta*sigma*sqrt(T))

		- pow(H/S,2.0*mu)*k_pr_CN(eta*y2 - eta*sigma*sqrt(T)));

	*Fcap = K*(pow(H/S,mu+lambda)*k_pr_CN(eta*z) 

		+ pow(H/S,mu - lambda)*k_pr_CN(eta*z - 2.0*eta*lambda*sigma*sqrt(T)));		

}	

double k_fmla_Barrier( 
					  //{	option type
					  //	0: call			1: put
					  //	0: up and effective	1: down and effective
					  //	0: in				1: out
					  long call_or_put, long up_or_down, long in_or_out,  
					  //}
					  //{	market parameter 
					  // S: spot
					  // X: strike
					  // H: barrier
					  // r: risk free interest rate
					  // div: continuous dividend rate
					  // sigma: volatility
					  // T: time to expiration
					  // K: rebate 	(KO case´Â at hit, KI case´Â at expiry)
					  double S, double X, double H, double r, double div, double sigma, double T, double K
					  //}
					  )
{

	//{	define constant  and variables	
	long CALL = 0;
	long PUT = 1;
	long UP = 0;
	long DOWN = 1;
	long INTYPE = 0;
	long OUTTYPE = 1;
	long alive_or_dead;
	long option_type;

	double eta, phi;
	double Acap, Bcap, Ccap, Dcap, Ecap, Fcap;		

	//}

	if ( T == 0.0 ) {	//	routines for T = 0

		if(in_or_out==INTYPE)	if(up_or_down==DOWN)	if(S<=H)	alive_or_dead = 1;
		else		alive_or_dead = 0;
		else					if(S>=H)	alive_or_dead = 1;
		else		alive_or_dead = 0;
		else				if(up_or_down==DOWN)	if(S<=H)	alive_or_dead = 0;
		else		alive_or_dead = 1;
		else					if(S>=H)	alive_or_dead = 0;
		else		alive_or_dead = 1;		

		if ( alive_or_dead == 0 ) return K;
		else {
			if ( call_or_put == 0 ) return max(0.0, S-X);
			else return max(0.0, X-S);
		}


	}

	//{	routines when spot is inside barrier range
	if( in_or_out == INTYPE && up_or_down == UP && call_or_put == CALL && S >= H ) 
		return k_fmla_BSCall(S, X, r, div, sigma, T);
	if( in_or_out == INTYPE && up_or_down == UP && call_or_put == PUT && S >= H ) 
		return k_fmla_BSPut(S, X, r, div, sigma, T);
	if( in_or_out == INTYPE && up_or_down == DOWN && call_or_put == CALL && S <= H ) 
		return k_fmla_BSCall(S, X, r, div, sigma, T);
	if( in_or_out == INTYPE && up_or_down == DOWN && call_or_put == PUT && S <= H ) 
		return k_fmla_BSPut(S, X, r, div, sigma, T);
	if( in_or_out == OUTTYPE && up_or_down == UP && S >= H ) return (K*exp(-r*T));
	if( in_or_out == OUTTYPE && up_or_down == DOWN && S <= H ) return (K*exp(-r*T)); 

	//}

	//{	indexing the option type
	if ( in_or_out == INTYPE )	if ( call_or_put == CALL )	if ( up_or_down == DOWN )	option_type = 1;
	else						option_type = 2;
	else						if ( up_or_down == DOWN )	option_type = 3;
	else						option_type = 4;
	else					if ( call_or_put == CALL )	if ( up_or_down == DOWN )	option_type = 5;
	else						option_type = 6;
	else						if ( up_or_down == DOWN )	option_type = 7;
	else						option_type = 8;
	//}

	switch ( option_type ) {	//	set eta and phi
	case 1:
		eta = 1.0;
		phi = 1.0;
		break;
	case 2:
		eta = -1.0;
		phi = 1.0;
		break;
	case 3:
		eta = 1.0;
		phi = -1.0;
		break;
	case 4:
		eta = -1.0;
		phi = -1.0;
		break;
	case 5:
		eta = 1.0;
		phi = 1.0;
		break;		
	case 6:
		eta = -1.0;
		phi = 1.0;
		break;
	case 7:
		eta = 1.0;
		phi = -1.0;
		break;
	case 8:
		eta = -1.0;
		phi = -1.0;
		break;
	}	

	//{	set intermediate variables
	set_intermediate_variables(	
		S,  X,  H,  r,  div,  sigma,  T,  K,		
		eta,  phi,																					
		&Acap,  &Bcap,  &Ccap,  &Dcap,  &Ecap,  &Fcap);	
	//}

	switch ( option_type ) {	//	pricing
	case 1:
		if ( X>H) return Ccap + Ecap;
		else return (Acap - Bcap + Dcap + Ecap);
		break;	
	case 2:
		if ( X>H) return Acap + Ecap;
		else return (Bcap - Ccap + Dcap + Ecap);	
		break;
	case 3:
		if ( X>H) return Bcap - Ccap + Dcap + Ecap;
		else return (Acap + Ecap);	
		break;
	case 4:
		if ( X>H) return Acap - Bcap + Dcap + Ecap;
		else return Ccap + Ecap;	
		break;
	case 5:
		if ( X>H) return Acap - Ccap + Fcap;
		else return Bcap - Dcap + Fcap;
		break;
	case 6:
		if ( X>H) return Fcap;
		else return (Acap - Bcap + Ccap - Dcap + Fcap);	
		break;
	case 7:
		if ( X>H) return Acap - Bcap + Ccap - Dcap + Fcap;
		else return Fcap;
		break;
	case 8:
		if ( X>H) return Bcap - Dcap + Fcap;
		else return Acap - Ccap + Fcap;
		break;
	}

	return 0.0;
}


double k_KI_CashOrNothing(double S, double H, double r, double div, double sigma, double T, double K, int UD)
{
	// Knock-In Cash (at expiration) or nothing option 
	// S: spot
	// H: barrier
	// r: risk free interest rate
	// div: continuous dividend rate
	// sigma: volatility
	// T: time to expiration
	// K: rebate 	
	// UD = 1: Up,  UD = 0: Down

	double b, mu, lambda;
	double x2, y2;
	double phi, eta;

	if (T == 0.0){
		if (UD == 1){
			if (S >= H) 
				return(K);
			else
				return(0.0);
		}
		else{
			if (S <= H)
				return(K);
			else
				return(0.0);
		}

	}


	if (UD == 1){
		if (S >= H) 
			return(K*exp(-r*T));
	}
	else{
		if (S <= H)
			return(K*exp(-r*T));
	}

	

	b = r - div;

	mu = (b - 0.5*sigma*sigma)/(sigma*sigma);

	lambda = sqrt(mu*mu + 2.0*r/(sigma*sigma));

	x2 = log(S/H)/(sigma*sqrt(T)) + (1.0 + mu)*sigma*sqrt(T);

	y2 = log(H/S)/(sigma*sqrt(T)) + (1.0 + mu)*sigma*sqrt(T);

	if (UD == 1){
		phi = 1.0;
		eta = -1.0;
	}
	else{
		phi = -1.0;
		eta = 1.0;
	}


	double B2, B4;

	B2 = K*exp(-r*T)*(k_pr_CN(phi*x2 - phi*sigma*sqrt(T)));


	B4 = K*exp(-r*T)*pow(H/S,2.0*mu)*k_pr_CN(eta*y2 - eta*sigma*sqrt(T));


	return(B2+B4);

}	

double k_KO_Rebate_at_Expiry(double S, double K, double H, double r, double div, double sigma, double T, double R, int UD, int CP)
{
	// KO option with rebate at expiry 

	double price;

	int up_and_down;
	int call_and_put;

	if (UD == 1) up_and_down = 0;
	else up_and_down = 1;

	if (CP == 1) call_and_put = 0;
	else call_and_put = 1;

	price = k_fmla_Barrier(call_and_put,up_and_down,1,S,K,H,r,div,sigma,T,0.0)
		+k_KI_CashOrNothing(S,H,r,div,sigma,T,R,UD);

	return(price);

}


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
					   )
{
	int n;
	double sum1, sum2;
	double d1,d2,d3,d4;
	double mu1,mu2,mu3;
	double b;
	double price;

	if ((S >= U) || (S <=L))
		return(0.0);

	if (T == 0.0){
		if (CP == 1)	
			return(max(S-K,0.0));
		else
			return(max(K-S,0.0));
	}

	b = r-q;
	sum1 = sum2 =0.0;
	if (CP == 1){
		for(n=-5; n <=5; n++){
			d1 = (log(S*k_power(U,2.0*n)/(K*k_power(L,2.0*n))) + (b+vol*vol*0.5)*T)/(vol*sqrt(T));
			d2 = (log(S*k_power(U,2.0*n)/(U*k_power(L,2.0*n))) + (b+vol*vol*0.5)*T)/(vol*sqrt(T));
			d3 = (log(k_power(L,2.0*n+2.0)/(K*S*k_power(U,2.0*n))) + (b+vol*vol*0.5)*T)/(vol*sqrt(T));
			d4 = (log(k_power(L,2.0*n+2.0)/(U*S*k_power(U,2.0*n))) + (b+vol*vol*0.5)*T)/(vol*sqrt(T));
			mu1 = 2.0*b/vol/vol+1.0;
			mu2 = 0.0;
			mu3 = 2.0*b/vol/vol+1.0;
			sum1 += k_power(k_power(U/L,1.0*n),mu1)*(k_pr_CN(d1)-k_pr_CN(d2)) -  k_power(k_power(L,n+1.0)/k_power(U,1.0*n)/S,mu3)*(k_pr_CN(d3)-k_pr_CN(d4));
			sum2 += k_power(k_power(U/L,1.0*n),mu1-2.0)*(k_pr_CN(d1-vol*sqrt(T))-k_pr_CN(d2-vol*sqrt(T)))
				- k_power(k_power(L,n+1.0)/k_power(U,1.0*n)/S,mu3-2)*(k_pr_CN(d3-vol*sqrt(T))-k_pr_CN(d4-vol*sqrt(T)));
			price = S*exp( (b-r)*T)*sum1-K*exp(-r*T)*sum2;
		}
	}
	else{
		for(n=-5; n <= 5; n++){
			d1 = (log(S*k_power(U,2.0*n)/k_power(L,2.0*n+1.0)) + (b+vol*vol*0.5)*T)/(vol*sqrt(T));
			d2 = (log(S*k_power(U,2.0*n)/(K*k_power(L,2.0*n))) + (b+vol*vol*0.5)*T)/(vol*sqrt(T));
			d3 = (log(k_power(L,2.0*n+2.0)/(L*S*k_power(U,2.0*n))) + (b+vol*vol*0.5)*T)/(vol*sqrt(T));
			d4 = (log(k_power(L,2.0*n+2.0)/(K*S*k_power(U,2.0*n))) + (b+vol*vol*0.5)*T)/(vol*sqrt(T));
			mu1 = 2.0*b/vol/vol+1.0;
			mu2 = 0.0;
			mu3 = 2.0*b/vol/vol+1.0;
			sum1 += k_power(k_power(U/L,1.0*n),mu1)*(k_pr_CN(d1)-k_pr_CN(d2)) -  k_power(k_power(L,n+1.0)/k_power(U,1.0*n)/S,mu3)*(k_pr_CN(d3)-k_pr_CN(d4));
			sum2 += k_power(k_power(U/L,1.0*n),mu1-2.0)*(k_pr_CN(d1-vol*sqrt(T))-k_pr_CN(d2-vol*sqrt(T)))
				- k_power(k_power(L,n+1.0)/k_power(U,1.0*n)/S,mu3-2)*(k_pr_CN(d3-vol*sqrt(T))-k_pr_CN(d4-vol*sqrt(T)));
			price = K*exp( -r*T)*sum2-S*exp((b-r)*T)*sum1;
		}
	}

	return(price);
}


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
						 )
{

	// To guarantee KO-KI parity

	return(  k_fmla_BSCall( S,K,r,q,vol,T) - k_DoubleBarrierKO(S,K,U,L,T,r,q,vol,CP));
	      
}


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
								 )
{
	int n;
	double sum1; //, sum2;
	double d1,d2,d3; //,d4;
	//double mu1,mu2,mu3;
	double b;
	//double price;
	double added;

	added = R*exp(-r*T);

	if ((S >= U) || (S <=L))
		return(added);

	if (T == 0.0){
		if ((S >= U) || (S <=L))
			return(R);
		else
			return(0.0);
	}

	b = r-q;


	double z;
	double alpha;
	double beta;
	double im;

	z = log(U/L);
	alpha = -0.5*(2.0*b/vol/vol-1.0);
	beta = -0.25*(2.0*b/vol/vol-1.0)*(2.0*b/vol/vol-1.0)-2.0*r/vol/vol;

	sum1 = 0.0;

	im = 1.0;

	for(n=1; n <100; n++){

		im *= (-1.0);

		d1 = 2.0*PI*n*R/z/z;
		d2 = ( k_power(S/L,alpha) - im*k_power(S/U,alpha) ) / (alpha*alpha + n*n*PI*PI/z/z);
		d3 = sin( n*PI/z *log(S/L))*exp(-0.5*( n*n*PI*PI/z/z -beta  )*vol*vol*T);

		sum1 += d1*d2*d3;
	}

	;

	return(added-sum1);
}

// One Touch Double Barrier KO (with rebate R at maturity)
double k_OneTouchDoubleBarrierKO(
								 double S,            // Spot
								 double U,            // Upper Barrier
								 double L,            // Lower Barrier
								 double T,            // time to maturity
								 double r,            // risk-free rate
								 double q,            // dividend rate
								 double vol,          // volatility
								 double R             // Rebate
								 )
{
	
	if (T == 0.0) return(R - k_OneTouchDoubleBarrierKI(S,U,L,T,r,q,vol,R) );

	return( R*exp(-r*T) - k_OneTouchDoubleBarrierKI(S,U,L,T,r,q,vol,R) );
}


double integrand_for_spread_option(double z, double K, double x1, double x2, double x11, double x22, double x12,
				 double b1, double b2, double nu, int CP)
// integrand for spread option approximation
{
	double temp1, temp2 ,temp3;

	double eps;

	if (CP == 1) eps = 1.0;
	else eps = -1.0;

	temp1 = nu * log( sqrt(x11)*exp(b2*z)/ (K  + x2*x12*exp(b1*z)/ x1/sqrt(x22) )  );
	temp2 = nu * log(  x1*x12*exp(b2*z)/x2/sqrt(x11)/ ( K + sqrt(x22)*exp(b1*z) ) );
	temp3 = nu * log( x1*x1 * exp(b2*z) / sqrt(x11) / (K + x2*x2/sqrt(x22)*exp(b1*z)));

	return(  (x1*k_pr_CN(temp1) - eps*x2*k_pr_CN(temp2) - eps* K*k_pr_CN(temp3))
		* k_pr_Npdf(0.0, 1.0, z));


}


double adaptive_quad_for_spread_option(double a, double b, double fa, double fb, double TOL,
				   double K, double x1, double x2, double x11, double x22, double x12,
				   double b1, double b2, double nu, int *level,int CP)
// adaptive quadrature for approximating spread option price 
{

	double h;
	double alpha, beta;
	double mll, ml, m, mr, mrr;
	double xx[5],y[5];
	double fmll, fml, fm, fmr, fmrr;
	double i1, i2;
	double Q;
	double xi;
	int i;

	int lev;

	lev = (*level) + 1;

	h = (b-a)/2; m= (a+b)/2;

	alpha = sqrt(2.0/3.0); beta = 1.0/sqrt(5.0);
	mll = m-alpha*h; ml = m-beta*h; mr = m+ beta*h; mrr = m+alpha*h;
	xx[0] = mll, xx[1] = ml, xx[2] = m, xx[3] = mr; xx[4] = mrr;



	for(i=0; i < 5; i++){
		xi = xx[i] ;
		y[i] = integrand_for_spread_option(xi,K,x1,x2,x11,x22,x12,b1,b2,nu,CP);
	}

	fmll = y[0]; fml = y[1]; fm = y[2]; fmr = y[3]; fmrr = y[4];

	i2 = h/6.0*(fa+fb+5.0*(fml+fmr));
	i1 = (h/1470.0)*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);


	if ( (fabs(i1-i2) < TOL) || ( mll <= a) || (b <= mrr) || (lev == 20) ){
		Q = i1;

	}
	else{
		Q = adaptive_quad_for_spread_option(a,mll,fa,fmll,TOL,K,x1,x2,x11,x22,x12,b1,b2,nu, &lev, CP);
		Q += adaptive_quad_for_spread_option(mll,ml,fmll,fml,TOL,K,x1,x2,x11,x22,x12,b1,b2,nu, &lev, CP);
		Q += adaptive_quad_for_spread_option(ml,m,fml,fm,TOL,K,x1,x2,x11,x22,x12,b1,b2,nu, &lev, CP);
		Q += adaptive_quad_for_spread_option(m,mr,fm,fmr,TOL,K,x1,x2,x11,x22,x12,b1,b2,nu, &lev, CP);
		Q += adaptive_quad_for_spread_option(mr,mrr,fmr,fmrr,TOL,K,x1,x2,x11,x22,x12,b1,b2,nu, &lev, CP);
		Q += adaptive_quad_for_spread_option(mrr,b,fmrr,fb,TOL,K,x1,x2,x11,x22,x12,b1,b2,nu, &lev, CP);

	}

	return(Q);

}


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
						  )
// Spread Option Approximation based on MUREX spread option formula
{
	double x1,x2,x11,x22,x12,b1,b2,nu;


	x1 = F1;
	x2 = F2;
	x11 = x1*x1*exp(v1*v1*T);
	x12 = x1*x2*exp(rho*v1*v2*T);
	x22 = x2*x2*exp(v2*v2*T);


	b1 = sqrt( log(x22/x2/x2)  );
	b2 = log(x12/x1/x2)/b1;
	nu = 1.0/sqrt(log(x11/x1/x1) - b2*b2 ); 


	double a,b,fa,fb;

	a = -5.0;
	b = 5.0;

	fa = integrand_for_spread_option(a,K,x1,x2,x11,x22,x12,b1,b2,nu,CP);
	fb = integrand_for_spread_option(b,K,x1,x2,x11,x22,x12,b1,b2,nu,CP); 

	int level;

	level = 1;

	double opValue;

	opValue = Disfactor * adaptive_quad_for_spread_option(a, b, fa, fb, 1.0e-8, K,x1,x2,x11,x22,x12,b1,b2,nu, &level,CP);


	return( opValue );
}



double call_iv(double dSpot,
					 double dStrike,
					 double dRate,
					 double dDividend,
					 double dVol,
					 double dExpiry,
					 double dCallPrice)
{
	double vol=dVol;
	double call;
	//call=call_price(100,0.011,100,365,0.1);
	//call=call_price(spot,r,k,t,vol);
	//call=call_vega(264.2,0.025,282.5,28,0.2);
	 
	for(int i=0;i<35;i++)
	{
		call=k_fmla_BSCall(dSpot,dStrike,dRate,dDividend,vol,dExpiry);
		vol = vol - (call-dCallPrice)/k_fmla_BSVega(dSpot,dStrike,dRate, dDividend, vol,dExpiry);
	}
 
	return vol;
}

double put_iv(double dSpot,
					 double dStrike,
					 double dRate,
					 double dDividend,
					 double dVol,
					 double dExpiry,
					 double dPutPrice)
{
	double vol=dVol;
	double put;
	//call=call_price(100,0.011,100,365,0.1);
	//call=call_price(spot,r,k,t,vol);
	//call=call_vega(264.2,0.025,282.5,28,0.2);
	 
	for(int i=0;i<35;i++)
	{
		put=k_fmla_BSPut(dSpot,dStrike,dRate,dDividend,vol,dExpiry);
		vol = vol - (put-dPutPrice)/k_fmla_BSVega(dSpot,dStrike,dRate, dDividend, vol,dExpiry);
	}
 
	return vol;
}































