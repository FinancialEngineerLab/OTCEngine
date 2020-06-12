////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	Program Name			: Zero Curve Murex Replication
// 	Version					: 0.0.0.0
// 	Author					: Hangseob CHO
// 	Date					: 2009.01.06
// 	Modified by				: 
// 	Modified at				:
// 	Copyright				: KDB Quant team
// 	Description				: Zero Curve Murex Replication
// 	Related Doc. Name		: Quant Report 2008-2
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef	INCLUDE_FILE_kdb_intRateTermStructure_hpp
#define	INCLUDE_FILE_kdb_intRateTermStructure_hpp

#include <stdio.h>
#include <stdlib.h>
#include <kdb_optimization.hpp>
#include <kdb_interpolation.hpp>
#include <kdb_date.hpp>
#include <kdb_cZeroCurve.hpp>
#include <kdb_miscellaneous.hpp>
#include <kdb_swap_exact_schedule.hpp>

////////////////////////////////////////////////////////////////////////////////////////////////////////
//	
//	The object of the below class runs Newton Method 
//	to find annually compounding zero rate that matches a swap rate of a given tenor,
//	provided with discount factors up to m-th swap pay date.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

class k_cNM_ZCMR: public k_opt_NewtonMethod
{

public:
	k_cNM_ZCMR( 
		double	*init_t, 
		double	*init_DF, 
		int		init_m, 
		int		init_n, 
		double	init_swap_rate 
	);
	~k_cNM_ZCMR();

	double	swap_rate;	//	the target rate that should be matched
	double 	*t;			//	time schedule of swap payment 
						//	~[0]: year fraction from today to swap start day
						//	~[i]: year fraction from today to i-th payment day
	double	*DF;		//  ~[i]: discount factor to today from i-th payment day

	int		m;			//	discount factors up to m-th payment day are given
	int		n;			//	number of payment day of the target swap

	double f( double x );
};

class k_zeroCurve_KRW_STD : public k_cZeroCurve
{

public:
	k_zeroCurve_KRW_STD( 
		date 	today,
		double	ON_rate, 
		double 	TN_rate, 
		double 	CD_rate,
		int		numSwap,
		int		*swap_tenor, 	//	[0]~[numSwap-1]	year * 4 ( i.e. num payment )
		double	*swap_rate	
	);
	k_zeroCurve_KRW_STD( 
		int		init_numZCDays,
		double	*init_ACZR_daily
	);
	~k_zeroCurve_KRW_STD();

	double ACZR( double t );
	
	int		numZCDays;			//	array size of ACZR_daily
	double	*ACZR_daily;		//	~[n-1] : Annually Compounding Zero Rate from today to n days later
};

class k_cNM_ZCMR_EUR: public k_opt_NewtonMethod
{

public:
	k_cNM_ZCMR_EUR( 
		double	*init_t, 
		double	*init_t2,        ///
		double	*init_DF, 
		int		init_m, 
		int		init_n, 
		double	init_swap_rate,
		double  dfTN              ///
		);
	~k_cNM_ZCMR_EUR();

	double	swap_rate;	//	the target rate that should be matched
	double 	*t;			//	time schedule of swap payment 
	double  *t0; 
	//	~[0]: year fraction from today to swap start day
	//	~[i]: year fraction from today to i-th payment day
	double	*DF;		//  ~[i]: discount factor to today from i-th payment day
	double DFTN;

	int		m;			//	discount factors up to m-th payment day are given
	int		n;			//	number of payment day of the target swap

	double f( double x );

};

class k_zeroCurve_EUR_STD : public k_cZeroCurve
{

public:
	k_zeroCurve_EUR_STD( 
		date 	today,
		double	ON_rate, 
		double 	TN_rate, 
		double 	W_rate,
		int     MnumSwap,
		int		numSwap,
		int		*swap_tenor, 	//	[0]~[numSwap-1]	year * 4 ( i.e. num payment )
		double	*Mswap_rate,
		double	*swap_rate	
		);
	k_zeroCurve_EUR_STD( 
		int		init_numZCDays,
		double	*init_ACZR_daily
		);
	~k_zeroCurve_EUR_STD();

	double ACZR( double t );

	int		numZCDays;			//	array size of ACZR_daily
	double	*ACZR_daily;		//	~[n-1] : Annually Compounding Zero Rate from today to n days later
};

class k_cZeroCurve_node_store : public k_cZeroCurve
{
public:

	k_cZeroCurve_node_store();
	~k_cZeroCurve_node_store();

	double	ACZR( double t );

	int		num_node;
	double	*t_node;		//	point of time node, whose size is num_node
	double	*ACZR_node;		//	annualy compounding zero rate value at node

};

void k_simple_fill_cZeroCurve_node_store( 
	   k_cZeroCurve_node_store   *pZC,
	   int						init__num_node,
	   double				   *init__t_node,
	   double				   *init__ACZR_node	
	   );

//	added at 090603
void k_IRTS_zeroCurve_from_daily_DF(
		k_cZeroCurve_node_store		*pZC,
		int							init__num_node,
		double						*DF_daily
		);

void k_add_end_node_cZeroCurve_node_store( 
									k_cZeroCurve_node_store   *pZC,
									double	new__t_node,
									double	new__ACZR_node
									);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//	addition for KRW_FX boot strap

class k_cNM_KRW_FX: public k_opt_NewtonMethod
{

public:
	k_cNM_KRW_FX( 
		k_cZeroCurve_node_store  *init__p_ZC_interm,
		k_cZeroCurve			   *init__p_ZC_USD_STD,	
		double					init__swap_rate,
		date					init__date__eval,
		date					init__date__mat,
		double					init__DF_KRW_settle_spot,
		double					init__DF_USD_swap_start_spot,
		double					init__pfixed_rate
		);
	~k_cNM_KRW_FX();

	k_cZeroCurve_node_store  *p_ZC_interm;
	k_cZeroCurve			   *p_ZC_USD_STD;
	double					swap_rate;	//	the target rate that should be matched
	date					date__eval;
	date					date__mat;
	double					DF_KRW_settle_spot;		//	DF of KRW_FX  over KRW settle ( + 1 SEOUL OPEN DAY) TO FX SPOT ( +2 SEOUL_NEW_YORK OPEN DAYS )
	double					DF_USD_swap_start_spot;//	DF of USD_STD over swap start TO FX spot
	double					pfixed_rate;	//	fixed rate for 1st USD floating payment
	
	double f( double x );
};

void k_KRW_FX_boot_strap(
					   //	target
					   k_cZeroCurve_node_store *p_ZC_KRW_FX, 
					   //	input
					   k_cZeroCurve			  *p_ZC_USD_STD,	

					   date	date__eval,

					   //	mkt quote
					   double	*rate__mkt
);

class k_cNM_USD_STD: public k_opt_NewtonMethod
{

public:
	k_cNM_USD_STD( 
		k_cZeroCurve_node_store  *init__p_ZC_interm,
		double					init__swap_rate,
		date					init__date__eval,
		date					init__date__swap_start,
		date					init__date__mat,
		double					init__pfixed_rate
		);
	~k_cNM_USD_STD();

	k_cZeroCurve_node_store  *p_ZC_interm;
	double					swap_rate;	//	the target rate that should be matched
	date					date__eval;
	date					date__swap_start;
	date					date__mat;
	double					pfixed_rate;	//	fixed rate for 1st USD floating payment

	double f( double x );
};

void k_USD_STD_boot_strap(
						//	target
						k_cZeroCurve_node_store *p_ZC_USD_STD, 
						//	INPUT: eval date
						date	date__eval,
						//	INPU: mkt quote
						double	*rate__mkt
						);

#endif