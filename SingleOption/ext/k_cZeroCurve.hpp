#ifndef	INCLUDE_FILE_kdb_cZeroCurve_hpp
#define	INCLUDE_FILE_kdb_cZeroCurve_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <kdb_optimization.hpp>
//#include <kdb_interpolation.hpp>
//#include <kdb_date.hpp>

//{ util functions before classes

double DF_from_ACZR( double r,  double t );
double ACZR_from_DF( double DF, double t );
double CR_from_DF( double DF, double t );
//}

class k_cZeroCurve
{

public:
	k_cZeroCurve();
	~k_cZeroCurve();

	virtual double ACZR( double t );		//	annual compounding zero rate

	// for funcions( t1, t2), they works when t1 == 0
	double DF( double t );					//	discount factor
	double DF( double t1, double t2 );			
	double CR( double t );					//	continuous compounding
	double CR( double t1, double t2 );
	double OR( double t, double year_fraction );		//	one-time compounding (like FRA)
	double OR( double t1, double t2, double year_fraction );
};

#endif