#include <kdb_cZeroCurve.hpp>

//{ util functions before classes

//////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Calculate discount foctor from annually compounding zero rate
//
/////////////////////////////////////////////////////////////////////////////////////////////////

double DF_from_ACZR( double r, double t )
{
	return exp( -t*log(1.+r) );
	//	correct even when t = 0.
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Calculate annually compounding zero rate from discount factor
//
/////////////////////////////////////////////////////////////////////////////////////////////////

double ACZR_from_DF( double DF, double t )
{
	return exp( -1./t*log(DF) ) - 1.;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Calculate continuously compounding zero rate from discount factor
//
/////////////////////////////////////////////////////////////////////////////////////////////////

double CR_from_DF( double DF, double t )
{
	return -log(DF)/t;
}

//}

//{ k_cZeroCurve class implementation
k_cZeroCurve::k_cZeroCurve()
{
}

k_cZeroCurve::~k_cZeroCurve()
{
}

double k_cZeroCurve::ACZR( double t)
{
	printf("error!!!: k_cZeroCurve: dummy function called");
	return 0.0;
}


/////////////////////////////////////////////////////////////////////////////////////////
//	
//	member function:
//		calculate Discount Factor
//	
////////////////////////////////////////////////////////////////////////////////////////

double k_cZeroCurve::DF( double t )
{
	return DF_from_ACZR( ACZR(t), t );
}

double k_cZeroCurve::DF( double t1, double t2 )
{
	//	value_at_t1 = disFactor(t1,t2)*value_at_t2
	
	//	DF(t2) = DF(t1)DF(t1,t2)
	//	DF(t1,t2) = DF(t2)/DF(t1)
	
	return DF(t2)/DF(t1);
}

////////////////////////////////////////////////////////////////////////////////////////
//	
//	member function:
//		calculate Continuous Compounding Zero Rate
//	
////////////////////////////////////////////////////////////////////////////////////////

double k_cZeroCurve::CR( double t )
{
	return log(1.+ACZR(t));
}

double k_cZeroCurve::CR( double t1, double t2 )
{
	//	DF = exp(-r*t)
	//	lnDF = -r*t
	//	-(lnDF)/t = r
	//	r = -ln(DF(t1,t2))/(t2-t1)
	//	r = -ln(DF(t2)/DF(t1))/(t2-t1)
	//	r = ( -lnDF(t2) + lnDF(t1) )/(t2-t1)
	//	r = ( r2*t2 - r1*t1 )/(t2-t1)
	return (CR(t2)*t2-CR(t1)*t1)/(t2-t1);
}

////////////////////////////////////////////////////////////////////////////////////////
//	
//	member function:
//		calculate One-time Compounding Zero Rate
//	
////////////////////////////////////////////////////////////////////////////////////////

double k_cZeroCurve::OR( double t, double year_fraction )
{
	//	1/DF = (1 + OR*year_fraction)
	//	OR*year_fraction = (1/DF -1)
	
	//	compounded once

	return (1./DF(t) -1.)/year_fraction;
}

double k_cZeroCurve::OR( double t1, double t2, double year_fraction )
{
	//	Forward Rate Agreement between t1 and t2
	//	compounded once

	//	(1 + OR*year_fraction) = 1/DF(t1,t2)
	//	OR*year_fraction = (1/DF(t1,t2) - 1.);
	return ( 1./DF(t1,t2) -1.)/year_fraction;
}

//}

