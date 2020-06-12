//	debug: USD와 JPY에 대해서도 구할 수 있게 바꿔준다.
//	현재는 KRW에 대해서만 bootstrapping이 일어난다.


#include <kdb_intRateTermStructure.hpp>

//{ cNM_ZCMR class implementation

/////////////////////////////////////////////////////////////////////////////////////////////////
//
//	constructor
//
/////////////////////////////////////////////////////////////////////////////////////////////////

k_cNM_ZCMR::k_cNM_ZCMR( double *init_t, double *init_DF, int init_m, int init_n, double init_swap_rate )
{
	int	i;
	
	swap_rate = init_swap_rate;
	
	t 	= init_t;
	
	m 	= init_m;
	n 	= init_n;
	
	DF = new double [n+1];
	
	for ( i = 0; i <= m; i++ )
		DF[i] = init_DF[i];
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//
//	destructor
//
/////////////////////////////////////////////////////////////////////////////////////////////////

k_cNM_ZCMR::~k_cNM_ZCMR()
{
	if ( DF ) delete[] DF;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Target fuction:
//		Newton method will find such x that makes f zero.
//
/////////////////////////////////////////////////////////////////////////////////////////////////

double k_cNM_ZCMR::f( double x ) 
{
	//	x is annually compounding zero rate at t[n]
	
	int	i;
	double *YF  = new double [n+1];		// ~[i] = Yearly Fraction from t[i-1] to t[i]
										//	i = 1, .. , n
	for ( i = 1; i <= n; i++ )
		YF[i] = ( t[i] - t[i-1] );
	
	double ACZR_m = ACZR_from_DF( DF[m], t[m] );
	
	DF[n] = DF_from_ACZR( x, t[n] );
	
	for ( i = m+1; i < n ; i++ ) {
		double ACZR = k_intp_linear( t[m], t[n], ACZR_m, x, t[i] );
		DF[i] = DF_from_ACZR( ACZR, t[i] );
	}
	
	double FCBP; 	//	fixed coupon bond price at t[0]
	
	FCBP = 0.;
	for ( i = 1; i <= n; i++ )
		FCBP += YF[i]*swap_rate*DF[i];
	FCBP += DF[n];
	
	return FCBP - DF[0];
}

k_cNM_ZCMR_EUR::k_cNM_ZCMR_EUR( double *init_t, double *init_t0, double *init_DF, int init_m, int init_n, double init_swap_rate, double dfTN )
{
	int	i;

	swap_rate = init_swap_rate;

	t 	= init_t;
	t0  = init_t0 ;

	m 	= init_m;
	n 	= init_n;

	DF = new double [n+1];
	DFTN = dfTN ;

	for ( i = 0; i <= m; i++ )
		DF[i] = init_DF[i];
}

k_cNM_ZCMR_EUR::~k_cNM_ZCMR_EUR()
{
	if ( DF ) delete[] DF;
}

double k_cNM_ZCMR_EUR::f( double x )   
{
	//	x is annually compounding zero rate at t[n]

	int	i;
	double *YF  = new double [n+1];		// ~[i] = Yearly Fraction from t[i-1] to t[i]
	//	i = 1, .. , n
	for ( i = 1; i <= n; i++ )
		YF[i] =   t0[i] - t0[i-1] ;
	

	double ACZR_m = ACZR_from_DF( DF[m], t[m]*360./365 );

	DF[n] = DF_from_ACZR( x, t[n]*360./365 );

	for ( i = m+1; i < n ; i++ ) {
		double ACZR = k_intp_linear( t[m]*360./365, t[n]*360./365, ACZR_m, x, t[i]*360./365 );
		DF[i] = DF_from_ACZR( ACZR, t[i]*360./365 );
	}

	double FCBP; 	//	fixed coupon bond price at t[0]

	

	FCBP = 0.;
	for ( i = 1; i <= n; i++ )
	{
		FCBP += YF[i]*swap_rate*DF[i]/DFTN ;
	}	
	FCBP += DF[n]/DFTN ;

	
	return FCBP -1.0 ; //- DF[0];
}


//}

//{ zeroCurveMurexReplicate class implementation

//////////////////////////////////////////////////////////////////////////////////////////////////
//	
//	constructor: type 01
//	
//	create zero curve from market rates ( overnight rate, tomorrow next rate, CD rate swap rates )
//	interpolation on annually compounding zero rate
//	
//////////////////////////////////////////////////////////////////////////////////////////////////

k_zeroCurve_KRW_STD::k_zeroCurve_KRW_STD( 
	date 	today,
	double	ON_rate,		//	overnight rate
	double 	TN_rate,		//	tomorrow next rate
	double 	CD_rate,		//	CD rate
	int		numSwap,		//
	int		*swap_tenor, 	//	[0]~[numSwap-1]	year * 4 ( i.e. num payment )
	double	*swap_rate	
)
{
	//{ memory and preliminary /////////////////////////////////////////////////////
	int		i;
	int		nb_pay = 20*4;

	double	 ACZR_TN;								//	from TN end date to today
	double	*ACZR_qtr	= new double [nb_pay+1];	//	~[0]: from ON end date ot today
													//	~[1]: from 3 mth to today		

	double	 DF_TN;									//	same as in ACZR
	double	*DF_qtr		= new double [nb_pay+1];

	double	 t_TN;									//	same as in ACZR
	double  *t_qtr 		= new double [nb_pay+1];	

	int		 tD_TN;
	int		*tD_qtr		= new int	 [nb_pay+1];
	
	vector<date> swap_pay_date;
	//}
	
	//{ schedule setup /////////////////////////////////////////////////////////////
	date ON_end = k_dt_workday(today+1,  1, SEOUL);
	date TN_end = k_dt_workday(ON_end+1, 1, SEOUL);
	
	swap_pay_date = k_dt_scheduleGenerator( ON_end, nb_pay , 3, 0, SEOUL );
		//	~[0       ]: first pay date
		//	~[nb_pay-1]: last  pay date
	
	tD_qtr[0] = k_dt_numDays( today, ON_end, 0 );
	tD_TN	  = k_dt_numDays( today, TN_end, 0 );
	for ( i = 1; i <= nb_pay; i++ ) 
		tD_qtr[i] = k_dt_numDays( today, swap_pay_date[i-1], 0 );
	
	
	t_TN = tD_TN/365.;
	for ( i = 0; i <= nb_pay; i++ )
		t_qtr[i] = tD_qtr[i]/365.;
	//}
	
	//{ boot strapping initial setup ///////////////////////////////////////////////
	DF_qtr[0] =           1./(1. + t_qtr[0]                     *ON_rate );
	DF_TN	  = DF_qtr[0]*1./(1. + k_dt_yearFrac( ON_end, TN_end, 0 )*TN_rate );
	DF_qtr[1] = DF_qtr[0]*1./(1. + (        t_qtr[1]-t_qtr[0]  )*CD_rate );
	
	ACZR_qtr[0] = ACZR_from_DF( DF_qtr[0], t_qtr[0] );
	ACZR_TN     = ACZR_from_DF( DF_TN    , t_TN     );
	ACZR_qtr[1] = ACZR_from_DF( DF_qtr[1], t_qtr[1] ); 
	//}
	
	//{ boot strapping /////////////////////////////////////////////////////////////
	int	i_swap;
	for ( i_swap = 0; i_swap < numSwap; i_swap++ ) {
		
		int m, n;
		
		m = ( i_swap == 0 )? 1: swap_tenor[i_swap -1];
		n = swap_tenor[i_swap];
		
		k_cNM_ZCMR oNM( t_qtr, DF_qtr, m, n, swap_rate[i_swap] );
		
		k_opt_NM_spec oNM_spec;
		oNM_spec.x0 		= swap_rate[i_swap];
		oNM_spec.dx 		= 1./100./100./100./100.;
		oNM_spec.err_tol 	= 1./100./100./100./100.;
		oNM_spec.max_itr 	= 1000;
		
		ACZR_qtr[n] 	= oNM.find_zero( oNM_spec );
		
		for ( i = m+1; i < n; i++ ) 
			ACZR_qtr[i] = k_intp_linear( t_qtr[m], t_qtr[n], ACZR_qtr[m], ACZR_qtr[n], t_qtr[i] );
		
		for ( i = m+1; i <= n; i++ )
			DF_qtr[i]   = DF_from_ACZR( ACZR_qtr[i], t_qtr[i] );
	}	
	//}

	//{ daily DF generation ////////////////////////////////////////////////////////
	numZCDays = tD_qtr[nb_pay] ;	
	
	ACZR_daily = new double [numZCDays];		// ~[0]: from day_1 to day_0 ( today )
	
	int	iD;

	for ( iD = 1; iD <= tD_qtr[0]; iD++ )
		ACZR_daily[iD-1] = ACZR_qtr[0]; 
	
	for ( iD = tD_qtr[0]+1; iD <= tD_TN; iD++ )
		ACZR_daily[iD-1] = k_intp_linear( t_qtr[0], t_TN, ACZR_qtr[0], ACZR_TN, iD/365. );
	
	for ( iD = tD_TN+1; iD <= tD_qtr[1]; iD++ )
		ACZR_daily[iD-1] = k_intp_linear( t_TN, t_qtr[1], ACZR_TN, ACZR_qtr[1], iD/365. );
	
	for ( i = 1; i < nb_pay; i++ ) {
		for ( iD = tD_qtr[i]+1; iD <= tD_qtr[i+1]; iD++ ) {
			ACZR_daily[iD-1] = k_intp_linear( t_qtr[i], t_qtr[i+1], ACZR_qtr[i], ACZR_qtr[i+1], iD/365. );
		}
	}
	
	//}
	
	//{	delete memory	/////////////////////////////////////////////////////////////
	if ( ACZR_qtr	) 		delete[]	ACZR_qtr	; 
	if ( DF_qtr   	)		delete[] 	DF_qtr 		;
	if ( t_qtr		)		delete[]	t_qtr		;
	if ( tD_qtr		)		delete[]	tD_qtr		;
	//}
}

//////////////////////////////////////////////////////////////////////////////////////////
//	
//	constructor: type 02
//	
//	create zero curve from Murex annually compounding zero rate
//	interpolation on annually compounding zero rate
//	
//////////////////////////////////////////////////////////////////////////////////////////

k_zeroCurve_KRW_STD::k_zeroCurve_KRW_STD( 
	int		init_numZCDays,
	double	*init_ACZR_daily	
)
{
	numZCDays = init_numZCDays;
	ACZR_daily = new double [numZCDays];
	
	for ( int i = 0; i < numZCDays; i++ )
		ACZR_daily[i] = init_ACZR_daily[i];
}	

//////////////////////////////////////////////////////////////////////////////////////////
//	
//	destructor
//	
//	delete DF_daily 
//	
//////////////////////////////////////////////////////////////////////////////////////////

k_zeroCurve_KRW_STD::~k_zeroCurve_KRW_STD()
{	
	if ( ACZR_daily )		delete[]	ACZR_daily;
}

//////////////////////////////////////////////////////////////////////////////////////////
//	
//	member function: 
//		calculate Annually Compounding Zero Rate
//	
//	It is linearly interpolated from ACZR_daily
//	
//////////////////////////////////////////////////////////////////////////////////////////

double	k_zeroCurve_KRW_STD::ACZR( double t )
{	
	int	LI = (int)floor(t*365.-1);
	
	if ( LI < 0 ) 				return ACZR_daily[0];
	if ( LI +1 > numZCDays-1 ) 	return ACZR_daily[numZCDays-1];
	return k_intp_linear( (double)LI, LI + 1., ACZR_daily[LI], ACZR_daily[LI+1], t*365.-1 );
}

//}

k_zeroCurve_EUR_STD::k_zeroCurve_EUR_STD( 
	date 	today,
	double	ON_rate,		//	overnight rate
	double 	TN_rate,		//	tomorrow next rate
	double 	W_rate,		//	CD rate
	int     MnumSwap,
	int		numSwap,		//
	int		*swap_tenor, 	//	[0]~[numSwap-1]	year * 4 ( i.e. num payment )
	double	*Mswap_rate, 
	double	*swap_rate	
	)
{
	//{ memory and preliminary /////////////////////////////////////////////////////
	int		i;
	int		Mnb_pay = 12 ;                 // From 1M to 12M
	int		nb_pay = 30 ;                  // From 1Y to 30Y

	double	 ACZR_TN;								//	from TN end date to today
	double   ACZR_W;
	double	*ACZR_mth	= new double [Mnb_pay+1];	//	~[0]: from ON end date to today
	double	*ACZR_qtr	= new double [nb_pay+1];	//	~[0]: from ON end date to today
	

	double	 DF_TN;                                 //	same as in ACZR
	double   DF_W ;
	double	*DF_mth		= new double [Mnb_pay+1];
	double	*DF_qtr		= new double [nb_pay+1];

	double	 t_TN;									//	same as in ACZR  // A/360
	double   t_W ;
    double  *t_mth 		= new double [Mnb_pay+1];
	double  *t_qtr 		= new double [nb_pay+1];	
	double  *t0_qtr 		= new double [nb_pay+1];	

	double  Date_Con = 360./365 ;                   // From A/360 To A/365  

	
	int		 tD_TN;
	int      tD_W ;
    int     *tD_mth     =  new int   [Mnb_pay+1];
	int		*tD_qtr		=  new int	 [nb_pay+1];

	vector<date> Mswap_pay_date;
	vector<date> swap_pay_date;
	//}

	//{ schedule setup /////////////////////////////////////////////////////////////
	date ON_end = k_dt_workday(today+1,  1, TARGET);                 ///** TARGET NOT LONDON
	date TN_end = k_dt_workday(ON_end+1, 1, TARGET);
	date W_end = k_dt_workday(TN_end+6+1, 1, TARGET);                ///** CEHCK THE MARKET CONVENTION 1W

	Mswap_pay_date = k_dt_scheduleGenerator( TN_end, Mnb_pay , 1, 0, TARGET );
	swap_pay_date = k_dt_scheduleGenerator( TN_end, nb_pay , 12, 0, TARGET );
	//	~[0       ]: first pay date
	//	~[nb_pay-1]: last  pay date

	
///** CHECK THE PAYMENT-DAY
//  printf("%%  %4d-%2d-%2d \n", ON_end.year, ON_end.month, ON_end.day);
// 	printf("%%  %4d-%2d-%2d \n", TN_end.year, TN_end.month, TN_end.day);
// 	printf("%%  %4d-%2d-%2d \n", W_end.year,  W_end.month,  W_end.day);
	
// 	for(i=0; i<=Mnb_pay-1 ; i++)
// 	  printf("!! %2d %4d-%2d-%2d \n",i, Mswap_pay_date[i].year, Mswap_pay_date[i].month, Mswap_pay_date[i].day);
// 	
// 	for(i=0; i<=nb_pay-1 ; i++)
// 	  printf("** %2d %4d-%2d-%2d \n", i, swap_pay_date[i].year, swap_pay_date[i].month, swap_pay_date[i].day);
	

	tD_TN	  = k_dt_numDays( today, TN_end, 0 );
	tD_W	  = k_dt_numDays( today, W_end, 0 );

	t_TN = tD_TN/360.;   t_W = tD_W/360.;     ///** FRACTION A/360

    
	tD_qtr[0] = k_dt_numDays( today, TN_end, 1 );

	for ( i = 1; i <= nb_pay; i++ ) 
		tD_qtr[i] = k_dt_numDays( today, swap_pay_date[i-1], 1 );

	for ( i = 0; i <= nb_pay; i++ ) 
		t0_qtr[i] = tD_qtr[i]/360.;

// 	for(i=1; i<=nb_pay-1 ; i++)
//  	{
//  	 		printf("%2d T0 %5d %11.8f %3d \n", i, tD_qtr[i], t0_qtr[i], tD_qtr[i]-tD_qtr[i-1]);
//  	}	

	tD_qtr[0] = k_dt_numDays( today, ON_end, 0 );
	tD_mth[0] = tD_qtr[0] ;
	

	for ( i = 1; i <= Mnb_pay; i++ ) 
		tD_mth[i] = k_dt_numDays( today, Mswap_pay_date[i-1], 0 );

	for ( i = 1; i <= nb_pay; i++ ) 
		tD_qtr[i] = k_dt_numDays( today, swap_pay_date[i-1], 0 );
	
	for ( i = 0; i <= Mnb_pay; i++ ) 
		t_mth[i] = tD_mth[i]/360.;

	for ( i = 0; i <= nb_pay; i++ ) 
		t_qtr[i] = tD_qtr[i]/360.;

// 	for(i=1; i<=nb_pay-1 ; i++)
// 	{
// 		printf("%2d %4d-%2d-%2d TA %5d \n", i, swap_pay_date[i].year, swap_pay_date[i].month, swap_pay_date[i].day
// 			,tD_qtr[i] );
// 	}	
	
	//{ boot strapping initial setup ///////////////////////////////////////////////
	DF_qtr[0] =           1./(1. + t_qtr[0]                     *ON_rate );
	DF_TN	  = DF_qtr[0]*1./(1. + k_dt_yearFrac( ON_end, TN_end, 2 )*TN_rate );      ///** FRACTION A/360
	DF_W	  = DF_TN*1./(1. + k_dt_yearFrac( TN_end, W_end, 2 )*W_rate );

    ACZR_qtr[0] = ACZR_from_DF( DF_qtr[0], t_qtr[0]*Date_Con );
	ACZR_TN     = ACZR_from_DF( DF_TN    , t_TN*Date_Con     );
	ACZR_W     = ACZR_from_DF( DF_W    , t_W*Date_Con     );
    	
    DF_mth[0] = DF_qtr[0] ;
	ACZR_mth[0] = ACZR_qtr[0] ;

// 	printf("ON %4d-%2d-%2d : %10.8f %10.8f\n", ON_end.year, ON_end.month, ON_end.day, DF_qtr[0],  ACZR_qtr[0]);
// 	printf("TN %4d-%2d-%2d : %10.8f %10.8f\n", TN_end.year, TN_end.month, TN_end.day, DF_TN, ACZR_TN);
// 	printf("W %4d-%2d-%2d : %10.8f %10.8f\n", W_end.year, W_end.month, W_end.day, DF_W, ACZR_W );
// 	printf("\n") ;

    //{ boot strapping initial setup ///////// WITHIN 1 YEAR USING SIMPLE IMTEREST ///
    double temp_rate ;

    for( i=1 ; i <=Mnb_pay ; i++ )
	{
       if(i==Mnb_pay)
	       temp_rate = swap_rate[0] ;
	   else
           temp_rate = Mswap_rate[i-1] ; 
	   
	   DF_mth[i] = DF_TN*1./(1. + ( t_mth[i]-t_TN  ) *temp_rate );
	   ACZR_mth[i] = ACZR_from_DF( DF_mth[i], t_mth[i] *Date_Con ); 

// 	   printf("%2d %4d-%2d-%2d : %10.8f %10.8f\n", i, 
// 		       Mswap_pay_date[i-1].year, Mswap_pay_date[i-1].month, Mswap_pay_date[i-1].day, DF_mth[i],  ACZR_mth[i]);
	}

	//{ boot strapping /////////////////////////////////////////////////////////////
	int	i_swap;

	t_qtr[1] = t_mth[12] ; ///**
	DF_qtr[1] = DF_mth[12] ; ///**
	ACZR_qtr[1] = ACZR_mth[12] ; ///**   

	for ( i_swap = 1; i_swap < numSwap ; i_swap++ ) {

		int m, n;

 		m = swap_tenor[i_swap -1];
 		n = swap_tenor[i_swap];

//		printf("\n\n M N %2d %2d \n", m,n) ;

      
//        for(i=m; i<=n;i++)
// 		   printf("%4d %10.8f :: %10.8f \n", tD_qtr[i], t_qtr[i], DF_qtr[i]) ;

		k_cNM_ZCMR_EUR oNM( t_qtr, t0_qtr, DF_qtr, m, n, swap_rate[i_swap], DF_TN );                             ///**

		k_opt_NM_spec oNM_spec;
		oNM_spec.x0 		= swap_rate[i_swap];

		oNM_spec.dx 		= 1./100./100./100./100./100. ;
		oNM_spec.err_tol 	= 1./100./100./100./100./100. ;
		oNM_spec.max_itr 	= 1000;

		ACZR_qtr[n] 	= oNM.find_zero( oNM_spec );

//		    printf("%10.8f %10.8f %10.8f %10.8f \n",t_qtr[m], t_qtr[n], ACZR_qtr[m], ACZR_qtr[n]) ;

		for ( i = m+1; i < n; i++ ) 
			ACZR_qtr[i] = k_intp_linear( t_qtr[m]*Date_Con, t_qtr[n]*Date_Con, ACZR_qtr[m], ACZR_qtr[n], t_qtr[i]*Date_Con );

		for ( i = m+1; i <= n; i++ )
			DF_qtr[i]   = DF_from_ACZR( ACZR_qtr[i], t_qtr[i]*Date_Con );
   
// 		for ( i = m; i <= n; i++ )
// 			printf("rate %12.10f dsf %12.10f\n",ACZR_qtr[i], DF_qtr[i]) ;


	}	
	//}

// 	for(i=0; i<=Mnb_pay ; i++)
// 		printf("%2d %10.8f  %10.8f\n",i, DF_mth[i], ACZR_mth[i] );
// 	
// 	for(i=0; i<=nb_pay ; i++)
// 	    printf("%2d %10.8f  %10.8f\n", i, DF_qtr[i], ACZR_qtr[i] );
	

	//{ daily DF generation ////////////////////////////////////////////////////////
	numZCDays = tD_qtr[nb_pay] ;	
//	printf("num %d \n",numZCDays) ;
 
	ACZR_daily = new double [numZCDays];		// ~[0]: from day_1 to day_0 ( today )

	int	iD;

	for ( iD = 1; iD <= tD_qtr[0]; iD++ )
		ACZR_daily[iD] = ACZR_qtr[0]; 

	for ( iD = tD_qtr[0]+1; iD <= tD_TN; iD++ )
		ACZR_daily[iD] = k_intp_linear( t_qtr[0], t_TN, ACZR_qtr[0], ACZR_TN, iD/360. );

	for ( iD = tD_TN+1; iD <= tD_W; iD++ )
		ACZR_daily[iD] = k_intp_linear( t_TN, t_W, ACZR_TN, ACZR_W, iD/360. );

	for ( iD = tD_W+1; iD <= tD_mth[1]; iD++ )
		ACZR_daily[iD] = k_intp_linear( t_W, t_mth[1], ACZR_W, ACZR_mth[1], iD/360. );

	for( i=1 ; i<=Mnb_pay-1 ; i++)
	  for ( iD = tD_mth[i]+1 ; iD <= tD_mth[i+1]; iD++ )
		 ACZR_daily[iD] = k_intp_linear( t_mth[i], t_mth[i+1], ACZR_mth[i], ACZR_mth[i+1], iD/360. );

	for( i=1 ; i<=nb_pay-1 ; i++)
	{
		for ( iD = tD_qtr[i]+1 ; iD <= tD_qtr[i+1]; iD++ )
		  {
         	ACZR_daily[iD] = k_intp_linear( t_qtr[i], t_qtr[i+1], ACZR_qtr[i], ACZR_qtr[i+1], iD/360. );
// 			if(iD==10095)
//               printf("% 5d %12.9f\n",iD,ACZR_daily[iD]*100.) ; 
		  }

	}


	//}

	//{	delete memory	/////////////////////////////////////////////////////////////
	if ( ACZR_qtr	) 		delete[]	ACZR_qtr	; 
	if ( DF_qtr   	)		delete[] 	DF_qtr 		;
	if ( t_qtr		)		delete[]	t_qtr		;
	if ( tD_qtr		)		delete[]	tD_qtr		;
	//}
}

k_zeroCurve_EUR_STD::k_zeroCurve_EUR_STD( 
	int		init_numZCDays,
	double	*init_ACZR_daily	
	)
{
	numZCDays = init_numZCDays;
	ACZR_daily = new double [numZCDays];

	for ( int i = 0; i < numZCDays; i++ )
		ACZR_daily[i] = init_ACZR_daily[i];
}	

k_zeroCurve_EUR_STD::~k_zeroCurve_EUR_STD()
{	
	if ( ACZR_daily )		delete[]	ACZR_daily;
}

double	k_zeroCurve_EUR_STD::ACZR( double t )
{	
	int	LI = (int)floor(t*365.-1);

	if ( LI < 0 ) 				return ACZR_daily[0];
	if ( LI +1 > numZCDays-1 ) 	return ACZR_daily[numZCDays-1];
	return k_intp_linear( (double)LI, LI + 1., ACZR_daily[LI], ACZR_daily[LI+1], t*365.-1 );
}





//{		cZeroCurve_node_store class implementation

k_cZeroCurve_node_store::k_cZeroCurve_node_store()
{
	t_node		= NULL;
	ACZR_node	= NULL;
}

k_cZeroCurve_node_store::~k_cZeroCurve_node_store()
{
	if ( t_node		)	delete[] t_node;
	if ( ACZR_node	)	delete[] ACZR_node;
	t_node		= NULL;
	ACZR_node	= NULL;
}

double k_cZeroCurve_node_store::ACZR( double t )
//	ACZR is calculated by interpolation from nearest two nodes.
//	If t < min{t_node} or max{t_node} < t, then flat extension.
{
	int	LI;		//	lower index s.t. t_node[LI] <= t < t_node[LI+1];

	LI	= k_misc_binsearch_lower_index( t_node, num_node, t );

	if ( LI == -1 ) return ACZR_node[0];
	else if ( LI == num_node ) return ACZR_node[num_node-1];
	else
	{
		return k_intp_linear( t_node[LI], t_node[LI+1], ACZR_node[LI], ACZR_node[LI+1], t);
	}
}

//}

//{		functions that fill cZeroCurve_daily_store or cZeroCurve_node_store
void k_simple_fill_cZeroCurve_node_store( 
	   k_cZeroCurve_node_store   *pZC,
	   int						init__num_node,
	   double				   *init__t_node,
	   double				   *init__ACZR_node	
	   )
{
	pZC->num_node	= init__num_node;
	pZC->t_node		= new double [init__num_node];
	pZC->ACZR_node	= new double [init__num_node];

	int		i;
	
	for ( i = 0; i < init__num_node; i++ )
	{
		pZC->t_node   [i] = init__t_node   [i];
		pZC->ACZR_node[i] = init__ACZR_node[i];
	}
}

//	added at 090603
void k_IRTS_zeroCurve_from_daily_DF(
		k_cZeroCurve_node_store		*pZC,
		int							init__num_node,
		double						*DF_daily
		)
//	DF_daily[i] : DF from today to (i+1)days later
//		i : 0 ~ init__num_node-1
{
	pZC->num_node	= init__num_node;
	pZC->t_node		= new double [init__num_node];
	pZC->ACZR_node	= new double [init__num_node];

	int		i;
	
	for ( i = 0; i < init__num_node; i++ )
	{
		pZC->t_node   [i] = (i+1)/365.;
		pZC->ACZR_node[i] = ACZR_from_DF( DF_daily[i], (i+1)/365. );
	}
}


void k_add_end_node_cZeroCurve_node_store( 
									k_cZeroCurve_node_store   *pZC,
									double	new__t_node,
									double	new__ACZR_node
									)
{
	double	*t_node		= new double [pZC->num_node+1];
	double	*ACZR_node	= new double [pZC->num_node+1];

	for ( int i = 0; i < pZC->num_node; i++ )
	{
		t_node[i]		= pZC->t_node   [i];
		ACZR_node[i]	= pZC->ACZR_node[i];
	}

	t_node	 [pZC->num_node] = new__t_node;
	ACZR_node[pZC->num_node] = new__ACZR_node;

	delete[] (pZC->t_node	);
	delete[] (pZC->ACZR_node);

	(pZC->num_node)++;
	pZC->t_node			= t_node;
	pZC->ACZR_node		= ACZR_node;
}

//}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//	addition for KRW_FX boot strap

k_cNM_KRW_FX::k_cNM_KRW_FX( 
					   k_cZeroCurve_node_store  *init__p_ZC_interm,
					   k_cZeroCurve			   *init__p_ZC_USD_STD,	
					   double					init__swap_rate,
					   date					init__date__eval,
					   date					init__date__mat,
					   double					init__DF_KRW_settle_spot,
					   double					init__DF_USD_swap_start_spot,
					   double					init__pfixed_rate
)
{
	p_ZC_interm			= init__p_ZC_interm;
	p_ZC_USD_STD		= init__p_ZC_USD_STD;
	swap_rate			= init__swap_rate;
	date__eval			= init__date__eval;
	date__mat			= init__date__mat;
	DF_KRW_settle_spot	= init__DF_KRW_settle_spot;
	DF_USD_swap_start_spot  = init__DF_USD_swap_start_spot;
	pfixed_rate			= init__pfixed_rate;
}

k_cNM_KRW_FX::~k_cNM_KRW_FX()
{
}

double k_cNM_KRW_FX::f( double x )
{
	//	return CRS NPV at FX spot date
	
	p_ZC_interm->ACZR_node[ p_ZC_interm->num_node - 1 ] = x;

	//	KRW_leg ( notional = 1 KRW )
	double PV_KRW_leg = k_swap_KRWCRS_FIXED_LEG( date__eval, date__eval, date__mat, S, swap_rate, 
		0, 0, SEOUL_NEW_YORK, p_ZC_interm);
	PV_KRW_leg /= DF_KRW_settle_spot;
	PV_KRW_leg += 1.;		//	exclude initial notional exchange

	//	USD_leg ( notional = 1 USD )
	double PV_USD_leg = 1./ DF_USD_swap_start_spot;
		//	PV at swap start is 1 for floating bond

	return PV_USD_leg - PV_KRW_leg;

	//	PVs are at USD-KRW spot date 
	//	KRWCRS_FIXED_LEG returns PV at KRW settlement date ( date_eval + 1 SEOUL OPEN DATE )
	//	KRW notional =  spot*USD_notional
	//	KRW PV of USD_leg = spot*PV_USD_leg
	//	so NPV = spot*USD_notioinal*( PV_USD_leg - PV_KRW_leg)
	
}

void k_KRW_FX_boot_strap(
					   //	target
					   k_cZeroCurve_node_store *p_ZC_KRW_FX, 
					   //	input
					   k_cZeroCurve			 *p_ZC_USD_STD,	

					   date	date__eval,

					   //	mkt quote
					   double	*rate__mkt
					   )
{
	int		i;
	int		nb__month_deposit	= 5;
	int		nb__swap			= 14;

	double	rate__ON				= rate__mkt[0]; 
	double 	rate__TN				= rate__mkt[1];
	double 	rate__week				= rate__mkt[2];
	double	*rate__month_deposit	= new double [nb__month_deposit];
	for ( i = 0; i < nb__month_deposit; i++ )
		rate__month_deposit[i] = rate__mkt[3+i];
	double	*rate__swap				= new double [nb__swap];
	for ( i = 0; i < nb__swap; i++ )
		rate__swap[i] = rate__mkt[3+nb__month_deposit+i];


	//	maturities are month
	int		mat__month_deposit[]	= {1,2,3,6,12};
	int		mat__swap[]			= {18,24,36,48,60,72,84,96,108,120,132,144,180,240};

	date	date__ON_end, date__TN_end, date__week_end;
	date	*date__month_deposit_end = new date [nb__month_deposit];
	date	*date__swap_end			 = new date [nb__swap];

	double	DF__ON_end, DF__TN_end, DF__week_end;
	double	*DF__month_deposit_end	= new double [nb__month_deposit];
	double	*DF__swap_end			= new double [nb__swap];

	//{	Set up dates	//////////////////////////////////////////////////////////////////////////////

	//	dates
	date__ON_end		= k_dt_workday(date__eval  +1, 1, SEOUL_NEW_YORK);	//	convention: following
	date__TN_end		= k_dt_workday(date__ON_end+1, 1, SEOUL_NEW_YORK);	//	convention: following
	date__week_end		= k_dt_workday(date__TN_end+7, 1, SEOUL_NEW_YORK);	//	convention: following

	for ( i = 0; i < nb__month_deposit; i++ )
	{
		date__month_deposit_end[i] = 
			k_dt_workday( k_dt_month_adder(date__TN_end, mat__month_deposit[i]), 0, SEOUL_NEW_YORK );
	}

	for ( i = 0; i < nb__swap; i++ )
	{
		date__swap_end[i] =
			k_dt_workday( k_dt_month_adder(date__TN_end, mat__swap[i]), 0, SEOUL_NEW_YORK );
	}

	//	t_ZR 을 이쪽에서 몽땅 계산해서 cNM에 넘겨주자.
	//	date 도 date__로 시작하는 변수로 몽땅 바꿔버릴까?

	double	t_ZR__ON_end,		      t_ZR__TN_end;	 
	double	t_ZR__week_end;			  
	double *t_ZR__month_deposit_end = new double [nb__month_deposit];
	double *t_ZR__swap_end			= new double [nb__swap];

	//	from date to t_ZR
	//		zero rate convention is act/365
	t_ZR__ON_end		= k_dt_yearFrac( date__eval, date__ON_end,		0 );
	t_ZR__TN_end		= k_dt_yearFrac( date__eval, date__TN_end,		0 );
	t_ZR__week_end		= k_dt_yearFrac( date__eval, date__week_end,	0 );
	
	for ( i = 0; i < nb__month_deposit; i++)
	{	
		t_ZR__month_deposit_end[i]		= k_dt_yearFrac( date__eval, date__month_deposit_end[i], 0 );
	}
	for ( i = 0; i < nb__swap; i++ )
	{
		t_ZR__swap_end[i]				= k_dt_yearFrac( date__eval, date__swap_end[i],         0 );
	}
	//}//////////////////////////////////////////////////////////////////////////////////////////////

	//{	calculate DF up to month deposit /////////////////////////////////////////////////////////////

	double dt;
	dt = k_dt_yearFrac( date__eval, date__ON_end, 0 );		//	Act/365
	DF__ON_end = 1./(1. + rate__ON*dt);
	dt = k_dt_yearFrac( date__ON_end, date__TN_end, 0 );	//	Act/365
	DF__TN_end = DF__ON_end*1./(1. + rate__TN*dt);
	dt = k_dt_yearFrac( date__TN_end, date__week_end, 0 );	//	Act/365
	DF__week_end = DF__TN_end*1./(1. + rate__week*dt);

	for ( i = 0; i < nb__month_deposit; i++ )
	{
		dt = k_dt_yearFrac( date__TN_end, date__month_deposit_end[i], 0 );
		DF__month_deposit_end[i] = DF__TN_end*1./(1. + rate__month_deposit[i]*dt);
	}

	//}//////////////////////////////////////////////////////////////////////////////////////////////

	//{	calculate ACZR up to month deposit /////////////////////////////////////////////////////////////

	double	ACZR_ON_end,		      ACZR_TN_end;	 
	double	ACZR_week_end;
	double *ACZR_month_deposit_end  = new double [nb__month_deposit];
	double *ACZR_swap_end			= new double [nb__swap];

	ACZR_ON_end		= ACZR_from_DF( DF__ON_end,   t_ZR__ON_end  );

	ACZR_TN_end		= ACZR_from_DF( DF__TN_end,   t_ZR__TN_end  );
	ACZR_week_end	= ACZR_from_DF( DF__week_end, t_ZR__week_end);

	for ( i = 0; i < nb__month_deposit; i++ )
	{
		ACZR_month_deposit_end[i] 
		= ACZR_from_DF( DF__month_deposit_end[i], t_ZR__month_deposit_end[i]);
	}

	//}///////////////////////////////////////////////////////////////////////////////////////////////

	//	generate zero curve with information up to month deposit
	
	double	*t_node_for_interm = new double [3+nb__month_deposit];
	double	*ACZR_for_interm   = new double [3+nb__month_deposit];

	t_node_for_interm[0] = t_ZR__ON_end;
	t_node_for_interm[1] = t_ZR__TN_end;
	t_node_for_interm[2] = t_ZR__week_end;
	for ( i = 0; i < nb__month_deposit; i++ )
		t_node_for_interm[3+i] = t_ZR__month_deposit_end[i];

	ACZR_for_interm  [0] = ACZR_ON_end;
	ACZR_for_interm  [1] = ACZR_TN_end;
	ACZR_for_interm  [2] = ACZR_week_end;
	for ( i = 0; i < nb__month_deposit; i++ )
		ACZR_for_interm  [3+i] = ACZR_month_deposit_end [i]; 
	
	k_cZeroCurve_node_store *p_ZC_interm	= new k_cZeroCurve_node_store;
	k_simple_fill_cZeroCurve_node_store
		( p_ZC_interm, 3+nb__month_deposit, t_node_for_interm, ACZR_for_interm);

	//	calc. fixing rate for USD floating first payment
	double yearFraction_1st_pay = k_dt_yearFrac(date__TN_end, date__month_deposit_end[3], 2); // act/360
	double DF_1st_pay = p_ZC_USD_STD->DF(t_ZR__TN_end, t_ZR__month_deposit_end[3]);

	//	(1+r*yearFraction) = 1/DF
	//	r = (1/DF -1)/yearFraction

	double pfixed_rate = (1./DF_1st_pay -1.)/yearFraction_1st_pay;

	date	date__FX_spot		= date__TN_end;
	date	date__KRW_settle	= k_dt_workday(date__eval  +1, 1, SEOUL);
	double	t_ZR__FX_spot		= (date__FX_spot    -date__eval)/365.;
	double  t_ZR__KRW_settle	= (date__KRW_settle -date__eval)/365.;

	double DF_KRW_settle_spot, DF_USD_swap_start_spot;

	DF_KRW_settle_spot = p_ZC_interm ->DF( t_ZR__KRW_settle, t_ZR__FX_spot );

	DF_USD_swap_start_spot = p_ZC_USD_STD->DF( t_ZR__TN_end, t_ZR__FX_spot );

	//	boot strap
	for ( i = 0; i < nb__swap; i++ )
	{
		//	add new node
		k_add_end_node_cZeroCurve_node_store( p_ZC_interm, t_ZR__swap_end[i], rate__swap[i] );

		k_cNM_KRW_FX	oNM_KRW_FX( p_ZC_interm, p_ZC_USD_STD, rate__swap[i], date__eval, 
								date__swap_end[i], DF_KRW_settle_spot, DF_USD_swap_start_spot,
								pfixed_rate													);

		k_opt_NM_spec oNM_spec;
		oNM_spec.x0 		= rate__swap[i];
		oNM_spec.dx 		= 1./100./100./100./100.;
		oNM_spec.err_tol 	= 1./100./100./100./100./100./100.;
		oNM_spec.max_itr 	= 1000;

		p_ZC_interm->ACZR_node[3+nb__month_deposit+i]	= oNM_KRW_FX.find_zero( oNM_spec );
	}

	k_simple_fill_cZeroCurve_node_store( p_ZC_KRW_FX, p_ZC_interm->num_node,
		p_ZC_interm->t_node, p_ZC_interm->ACZR_node);

	delete p_ZC_interm;

	//	delete
	delete[]	rate__month_deposit;
	delete[]	rate__swap;
	delete[]	date__month_deposit_end;
	delete[]	date__swap_end;
	delete[]	DF__month_deposit_end;
	delete[]	DF__swap_end;
	delete[]	t_ZR__month_deposit_end;
	delete[]	t_ZR__swap_end;
	delete[]	ACZR_month_deposit_end;
	delete[]	ACZR_swap_end;
	delete[]	t_node_for_interm;
	delete[]	ACZR_for_interm;
}

k_cNM_USD_STD::k_cNM_USD_STD( 
						 k_cZeroCurve_node_store  *init__p_ZC_interm,
						 double					init__swap_rate,
						 date					init__date__eval,
						 date					init__date__swap_start,
						 date					init__date__mat,
						 double					init__pfixed_rate
						 )
{
	p_ZC_interm			= init__p_ZC_interm;
	swap_rate			= init__swap_rate;
	date__eval			= init__date__eval;
	date__swap_start	= init__date__swap_start;
	date__mat			= init__date__mat;
	pfixed_rate			= init__pfixed_rate;
}

k_cNM_USD_STD::~k_cNM_USD_STD()
{
}

double k_cNM_USD_STD::f( double x )
{
	//	return USD IRS NPV 

	p_ZC_interm->ACZR_node[ p_ZC_interm->num_node - 1 ] = x;

	// PV is PV at swap start date 

	//	PV fixed coupon 
	double PV_fixed_coupon_bond 
		= k_swap_USDIRS_FIXED_LEG( date__eval, date__eval, date__mat, A, swap_rate, 0, 2, p_ZC_interm);

	//	PV fixed coupon bond
	double	t_ZR__swap_start = (date__swap_start - date__eval)/365.;
	double	t_ZR__mat		 = (date__mat        - date__eval)/365.;

	PV_fixed_coupon_bond += p_ZC_interm->DF(t_ZR__swap_start, t_ZR__mat);

	return PV_fixed_coupon_bond-1;
}

void k_USD_STD_boot_strap(
						//	target
						k_cZeroCurve_node_store *p_ZC_USD_STD, 
						//	INPUT: eval date
						date	date__eval,
						//	INPU: mkt quote
						double	*rate__mkt
						)
{
	//	read mkt_rate
	double	rate__ON	= rate__mkt[0]; 
	double 	rate__TN	= rate__mkt[1]; 
	double 	rate__week	= rate__mkt[2];

	double	rate__month_deposit[12];
	double	rate__swap[14];

	int		nb__month_deposit	= 12;
	int		nb__swap			= 14;

	int		i;

	for( i = 0; i < nb__month_deposit; i++ )
	{
		rate__month_deposit[i] = rate__mkt[3+i];
	}

	for( i = 0; i < nb__swap; i++ )
	{
		rate__swap[i] = rate__mkt[3+nb__month_deposit+i];
	}

	//	maturities are month
	int		mat__month_deposit[]	= {1,2,3,4,5,6,7,8,9,10,11,12};
	int		mat__swap[]			= {24,36,48,60,72,84,96,108,120,144,180,240,300,360};

	date	date__ON_end, date__TN_end, date__week_end;
	date	*date__month_deposit_end = new date [nb__month_deposit];
	date	*date__swap_end			 = new date [nb__swap];

	double	DF__ON_end, DF__TN_end, DF__week_end;
	double	*DF__month_deposit_end	= new double [nb__month_deposit];
	double	*DF__swap_end			= new double [nb__swap];

	//{	Set up dates	//////////////////////////////////////////////////////////////////////////////

	//	dates
	date__ON_end		= k_dt_workday(date__eval  +1, 1, NEW_YORK);	//	convention: following
	date__TN_end		= k_dt_workday(date__ON_end+1, 1, NEW_YORK);	//	convention: following
	date__week_end		= k_dt_workday(date__TN_end+7, 1, NEW_YORK);	//	convention: following

	date date__swap_start 
		= k_dt_workday( k_dt_workday(date__eval+1, 1, NEW_YORK_LONDON)+1, 1, NEW_YORK_LONDON); 

	for ( i = 0; i < nb__month_deposit; i++ )
	{
		date__month_deposit_end[i] = 
			k_dt_workday( k_dt_month_adder(date__TN_end, mat__month_deposit[i]), 0, NEW_YORK );
	}

	for ( i = 0; i < nb__swap; i++ )
	{
		date__swap_end[i] =
			k_dt_workday( k_dt_month_adder(date__swap_start, mat__swap[i]), 0, NEW_YORK_LONDON );
	}

	//	t_ZR 을 이쪽에서 몽땅 계산해서 cNM에 넘겨주자.
	//	date 도 date__로 시작하는 변수로 몽땅 바꿔버릴까?

	double	t_ZR__ON_end,		      t_ZR__TN_end;	 
	double	t_ZR__week_end;			  
	double *t_ZR__month_deposit_end = new double [nb__month_deposit];
	double *t_ZR__swap_end			= new double [nb__swap];

	//	from date to t_ZR
	//		zero rate convention is act/365
	t_ZR__ON_end		= k_dt_yearFrac( date__eval, date__ON_end,		0 );
	t_ZR__TN_end		= k_dt_yearFrac( date__eval, date__TN_end,		0 );
	t_ZR__week_end		= k_dt_yearFrac( date__eval, date__week_end,	0 );

	for ( i = 0; i < nb__month_deposit; i++)
	{	
		t_ZR__month_deposit_end[i]		= k_dt_yearFrac( date__eval, date__month_deposit_end[i], 0 );
	}
	for ( i = 0; i < nb__swap; i++ )
	{
		t_ZR__swap_end[i]				= k_dt_yearFrac( date__eval, date__swap_end[i],         0 );
	}
	//}//////////////////////////////////////////////////////////////////////////////////////////////

	//{	calculate DF up to month deposit /////////////////////////////////////////////////////////////

	double dt;
	dt = k_dt_yearFrac( date__eval, date__ON_end, 2 );		//	Act/360
	DF__ON_end = 1./(1. + rate__ON*dt);
	dt = k_dt_yearFrac( date__ON_end, date__TN_end, 2 );	//	Act/360
	DF__TN_end = DF__ON_end*1./(1. + rate__TN*dt);
	dt = k_dt_yearFrac( date__TN_end, date__week_end, 2 );	//	Act/360
	DF__week_end = DF__TN_end*1./(1. + rate__week*dt);

	for ( i = 0; i < nb__month_deposit; i++ )
	{
		dt = k_dt_yearFrac( date__TN_end, date__month_deposit_end[i], 2 );
		DF__month_deposit_end[i] = DF__TN_end*1./(1. + rate__month_deposit[i]*dt);
	}

	//}//////////////////////////////////////////////////////////////////////////////////////////////

	//{	calculate ACZR up to month deposit /////////////////////////////////////////////////////////////

	double	ACZR_ON_end,		      ACZR_TN_end;	 
	double	ACZR_week_end;
	double *ACZR_month_deposit_end  = new double [nb__month_deposit];
	double *ACZR_swap_end			= new double [nb__swap];

	ACZR_ON_end		= ACZR_from_DF( DF__ON_end,   t_ZR__ON_end  );

	ACZR_TN_end		= ACZR_from_DF( DF__TN_end,   t_ZR__TN_end  );
	ACZR_week_end	= ACZR_from_DF( DF__week_end, t_ZR__week_end);

	for ( i = 0; i < nb__month_deposit; i++ )
	{
		ACZR_month_deposit_end[i] 
		= ACZR_from_DF( DF__month_deposit_end[i], t_ZR__month_deposit_end[i]);
	}

	//}///////////////////////////////////////////////////////////////////////////////////////////////

	//	generate zero curve with information up to month deposit

	double	*t_node_for_interm = new double [3+nb__month_deposit];
	double	*ACZR_for_interm   = new double [3+nb__month_deposit];

	t_node_for_interm[0] = t_ZR__ON_end;
	t_node_for_interm[1] = t_ZR__TN_end;
	t_node_for_interm[2] = t_ZR__week_end;
	for ( i = 0; i < nb__month_deposit; i++ )
		t_node_for_interm[3+i] = t_ZR__month_deposit_end[i];

	ACZR_for_interm  [0] = ACZR_ON_end;
	ACZR_for_interm  [1] = ACZR_TN_end;
	ACZR_for_interm  [2] = ACZR_week_end;
	for ( i = 0; i < nb__month_deposit; i++ )
		ACZR_for_interm  [3+i] = ACZR_month_deposit_end [i]; 

	k_cZeroCurve_node_store *p_ZC_interm	= new k_cZeroCurve_node_store;
	k_simple_fill_cZeroCurve_node_store
		( p_ZC_interm, 3+nb__month_deposit, t_node_for_interm, ACZR_for_interm);

	//	calc. fixing rate for USD floating first payment
	date date__first_floating_end = 
		k_dt_workday( k_dt_month_adder(date__swap_start, 3), 0, NEW_YORK_LONDON );

	double yearFraction_1st_pay 
		= k_dt_yearFrac(date__swap_start, date__first_floating_end, 2); // act/360

	double t_ZR__swap_start, t_ZR__first_floating_end;
	t_ZR__swap_start			= (date__swap_start         - date__eval)/365.;
	t_ZR__first_floating_end	= (date__first_floating_end - date__eval)/365.;

	double DF_1st_pay = p_ZC_interm->DF(t_ZR__swap_start, t_ZR__first_floating_end);

	//	(1+r*yearFraction) = 1/DF
	//	r = (1/DF -1)/yearFraction

	double pfixed_rate = (1./DF_1st_pay -1.)/yearFraction_1st_pay;

	//	boot strap
	for ( i = 0; i < nb__swap; i++ )
	{
		//	add new node
		k_add_end_node_cZeroCurve_node_store( p_ZC_interm, t_ZR__swap_end[i], rate__swap[i] );

		k_cNM_USD_STD	oNM_USD_STD( p_ZC_interm, rate__swap[i], date__eval, date__swap_start, 
			date__swap_end[i], pfixed_rate								);

		k_opt_NM_spec oNM_spec;
		oNM_spec.x0 		= rate__swap[i];
		oNM_spec.dx 		= 1./100./100./100./100.;
		oNM_spec.err_tol 	= 1./100./100./100./100./100./100.;
		oNM_spec.max_itr 	= 1000;

		p_ZC_interm->ACZR_node[3+12+i] = oNM_USD_STD.find_zero( oNM_spec );

	}

	k_simple_fill_cZeroCurve_node_store( p_ZC_USD_STD, p_ZC_interm->num_node,
		p_ZC_interm->t_node, p_ZC_interm->ACZR_node);

	delete p_ZC_interm;

	//	delete
	delete[]	date__month_deposit_end;
	delete[]	date__swap_end;
	delete[]	DF__month_deposit_end;
	delete[]	DF__swap_end;
	delete[]	t_ZR__month_deposit_end;
	delete[]	t_ZR__swap_end;
	delete[]	ACZR_month_deposit_end;
	delete[]	ACZR_swap_end;
	delete[]	t_node_for_interm;
	delete[]	ACZR_for_interm;
}
