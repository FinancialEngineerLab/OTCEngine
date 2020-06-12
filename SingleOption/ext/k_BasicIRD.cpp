#include "kdb_BasicIRD.hpp"
// #inlcude "kdb_date.hpp" && "kdb_swap_exact_schedule.hpp"


double k_Black_d1( double F_0, double K, double T, double sigma )
// parameter for Black-Scholes Formula 
{	
	//	F_0	:	forward price (maturity T, seen at t=0)
	//	K	:	strike
	//	T	:	option maturity 
	//	sigma	:	volatility of underlying
	return (log(F_0/K) + sigma*sigma*T/2.0)/(sigma*sqrt(T));
}	

double k_BlackCall( double disBond_0_Tstar, double F_0, double K, double T, double sigma)
// Black-Scholes Call Formula 
{
	//	Tstar	:	payoff delivery time
	//	disBond_0_Tstar	:	zero-coupon bond price maturing at Tstar
	//	F_0	:	forward price (maturity T, seen at t=0)
	//	K	:	strike
	//	T	:	option maturity
	//	sigma	:	volatility of underlying 

	double d1 = k_Black_d1( F_0, K, T, sigma );
	double d2 = d1 - sigma*sqrt(T);

	return disBond_0_Tstar*(F_0*k_pr_CN(d1) - K*k_pr_CN(d2));
}

double k_BlackPut( double disBond_0_Tstar, double F_0, double K, double T, double sigma)
// Black-Scholes Put Formula 
{
	//	Tstar	:	payoff delivery time
	//	disBond_0_Tstar	:	zero-coupon bond price maturing at Tstar
	//	F_0	:	forward price (maturity T, seen at t=0)
	//	K	:	strike
	//	T	:	option maturity
	//	sigma	:	volatility of underlying

	double d1 = k_Black_d1( F_0, K, T, sigma );
	double d2 = d1 - sigma*sqrt(T);

	return disBond_0_Tstar*(K*k_pr_CN(-d2) - F_0*k_pr_CN(-d1));
}

double k_caplet( double disBond_mat, double F, double K, double c_s, double c_int, double L, double sigma)
{
	//	disBond_mat	:	disFactor from effective date to payment date                  
	//                     or Zero Coupon Bond price from effective date to payment date
	//  F           :   today's forward rate 
	//	K			:	strike
	//	c_s			:	fixing day                        
	//	c_int	    :	calculation period size = calculation end - calculation start 
	//											= payment day - caplet maturity
	//                                            considering day-count convention
	//	L			:	principal
	//	sigma		:	volatility
	return L*c_int*k_BlackCall(disBond_mat, F, K, c_s, sigma);
}

double k_floorlet( double disBond_mat, double F, double K, double c_s, double c_int, double L, double sigma)
{
	//	disBond_mat	:	disFactor from effective date to payment date                  
	//                     or Zero Coupon Bond price from effective date to payment date
	//  F           :   today's forward rate 
	//	K			:	strike
	//	c_s			:	fixing day                        
	//	c_int	    :	calculation period size = calculation end - calculation start 
	//											= payment day - floorlet maturity
	//                                            considering day-count convention
	//	L			:	principal
	//	sigma		:	volatility
	return L*c_int*k_BlackPut(disBond_mat, F, K, c_s, sigma);
}



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
			 )
{
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


	vector<date> payment;
	int *pdate;
	int *fdate;
	date *fixing;
	int num_rol;
	int i;
	int num_month;
		
	
	if (ptype == M)
		num_month = 1;
	else if (ptype == Q)
		num_month = 3;
	else if (ptype == S)
		num_month = 6;
	else if (ptype == A)
		num_month = 12;
	else 
		return(-2000.0); /// payment type error

	num_rol = ( (capmat.year-today.year)*12 + (capmat.month-today.month))/num_month ;

    date effective_date;

	int edate;

	effective_date = k_swap_getSpotDate(today,pay_cal);
	edate = k_dt_date_to_count(effective_date) - k_dt_date_to_count(today);
    payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,pay_cal);
	
	pdate = (int *)malloc(sizeof(int)*num_rol);   // payment date = calculation end 
	fdate = (int *)malloc(sizeof(int)*num_rol);   // fixing date 
	fixing = (date *)malloc(sizeof(date)*num_rol);
   
	for(i=0; i < num_rol; i++){
		pdate[i] = k_dt_date_to_count(payment[i]) - k_dt_date_to_count(today);
		fixing[i] = k_swap_fix_date_shift(payment[i],fix_cal);
		fdate[i] = k_dt_date_to_count(fixing[i]) - k_dt_date_to_count(today);
	}
   	

	double price;
	int frdate;

	price = 0.0;

	vector<date> fraed;
	date fstart;
	int fsdate;

	for (i= 1; i< num_rol; i++) {
	
		fstart = k_dt_getSpotDate(fixing[i-1],pay_cal);
		fraed = k_dt_scheduleGenerator(fstart,1,num_month,BizDayConv,pay_cal);

		fsdate = k_dt_date_to_count(fstart) - k_dt_date_to_count(today);
		frdate = k_dt_date_to_count(fraed[0]) - k_dt_date_to_count(today);

		price += k_caplet( (*zCurve).DF(edate/365.0,pdate[i]/365.0),
			(*zCurve).OR(fsdate/365.0,frdate/365.0, k_dt_yearFrac(fstart,fraed[0],DayCountConv)),
			strike, fdate[i-1]/365.0, k_dt_yearFrac(payment[i-1],payment[i],DayCountConv),L,sigma);

         
		//printf("caplet %d %.5f %.5f %.5f %.5f \n,")

	//	printf("%d %.10f %d %d %d %.15f %.10f\n",i,100.0*(*zCurve).OR(fsdate/365.0,frdate/365.0, k_dt_yearFrac(fstart,fraed[0],DayCountConv)),
		//	k_dt_numDays(fstart,fraed[0],DayCountConv),k_dt_date_to_count(fstart),
		//	k_dt_date_to_count(fraed[0]),(*zCurve).DF(fsdate/365.0),(*zCurve).DF(edate/365.0,pdate[i]/365.0));
	
	}

	if (pdate) free(pdate);
    if (fdate) free(fdate);
	if (fixing) free(fixing);
    
	return price;	
}


double k_cap( k_cZeroCurve *zCurve, double strike, int nb_fixing, double L, double sigma, PayType ptype )
{
	 // simplified cap pricer 
	 //	return price of cap using flat vol
	 //	the caplet whose payoff fixed at today is not included

 	 //	nb_fixing : number of fixings = number of caplets, 	
	 //	L : principal
	 //	sigma: volatility
	 //  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually

     int k;

	 double price = 0.0;
	 double tau;

	 if (ptype == M)
		 tau = 1.0/12;
	 else if (ptype == Q)
		 tau = 0.25;
	 else if (ptype == S)
		 tau = 0.5;
	 else if (ptype == A)
		 tau = 1.0;

	 for (k= 1; k<= nb_fixing; k++) {

		 price +=  k_caplet( (*zCurve).DF(tau*(k+1)), (*zCurve).OR(tau*k,tau*(k+1),tau), strike, tau*k, tau, L, sigma);   
		         
	 }
	 return price;	
}


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
			 )
{
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


	vector<date> payment;
	int *pdate;
	int *fdate;
	date *fixing;
	int num_rol;
	int i;
	int num_month;


	if (ptype == M)
	  num_month = 1;
	else if (ptype == Q)
	 num_month = 3;
	else if (ptype == S)
	 num_month = 6;
	else if (ptype == A)
	 num_month = 12;
	else 
	 return(-2000.0); /// payment type error

	num_rol = ( (floormat.year-today.year)*12 + (floormat.month-today.month))/num_month ;

    date effective_date;

	int edate;
	effective_date = k_swap_getSpotDate(today,pay_cal);
	edate = k_dt_date_to_count(effective_date) - k_dt_date_to_count(today);
	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,pay_cal);

	pdate = (int *)malloc(sizeof(int)*num_rol);   // payment date = calculation end 
	fdate = (int *)malloc(sizeof(int)*num_rol);   // fixing date 
	fixing = (date *)malloc(sizeof(date)*num_rol);

	for(i=0; i < num_rol; i++){
	 pdate[i] = k_dt_date_to_count(payment[i]) - k_dt_date_to_count(today);
	 fixing[i] = k_swap_fix_date_shift(payment[i],fix_cal);
	 fdate[i] = k_dt_date_to_count(fixing[i]) - k_dt_date_to_count(today);
	}


	double price;
	int frdate;
	price = 0.0;

	vector<date> fraed;
	date fstart;
	int fsdate;

	for (i= 1; i< num_rol; i++) {

	 fstart = k_dt_getSpotDate(fixing[i-1],pay_cal);
	 fraed = k_dt_scheduleGenerator(fstart,1,num_month,BizDayConv,pay_cal);
	 fsdate = k_dt_date_to_count(fstart) - k_dt_date_to_count(today);
	 frdate = k_dt_date_to_count(fraed[0]) - k_dt_date_to_count(today);

	 price += k_floorlet( (*zCurve).DF(edate/365.0,pdate[i]/365.0),
			 (*zCurve).OR(fsdate/365.0,frdate/365.0, k_dt_yearFrac(fstart,fraed[0],DayCountConv)),
 		     strike, fdate[i-1]/365.0, k_dt_yearFrac(payment[i-1],payment[i],DayCountConv),L,sigma);

	// printf("%d %.10f %d %d %d %.15f %.10f\n",i,100.0*(*zCurve).OR(fsdate/365.0,frdate/365.0, k_dt_yearFrac(fstart,fraed[0],DayCountConv)),
		//	 k_dt_numDays(fstart,fraed[0],DayCountConv),k_dt_date_to_count(fstart),
		//	 k_dt_date_to_count(fraed[0]),(*zCurve).DF(fsdate/365.0),(*zCurve).DF(edate/365.0,pdate[i]/365.0));

	 }

	if (pdate) free(pdate);
    if (fdate) free(fdate);
    if (fixing) free(fixing);

	 return price;	
}


double k_floor( k_cZeroCurve *zCurve, double strike, int nb_fixing, double L, double sigma, PayType ptype )
{
	 // simplified floor pricer 
	 //	return price of floor using flat vol
	 //	the floorlet whose payoff fixed at today is not included

     //	nb_fixing : number of fixings = number of floorlets, 	
     //	L : principal
     //	sigma: volatility
	 //  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually

	 int k;

	 double price = 0.0;
	 double tau;

	 if (ptype == M)
		 tau = 1.0/12;
	 else if (ptype == Q)
		 tau = 0.25;
	 else if (ptype == S)
		 tau = 0.5;
	 else if (ptype == A)
		 tau = 1.0;

	 for (k= 1; k<= nb_fixing; k++) {
			 price +=  k_floorlet( (*zCurve).DF(tau*(k+1)), (*zCurve).OR(tau*k,tau*(k+1.0),tau), strike, tau*k, tau, L, sigma);   
	 }
	 return price;	
}

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
						 ){
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


    if ( k_dt_date_to_count(today) > k_dt_date_to_count(omat))
		return(-1000.0);   // cannot evaluate forward swap rate

//	if ( k_dt_date_to_count(today) == k_dt_date_to_count(omat))
//		return( (*zCurve).ACZR( (k_dt_date_to_count(smat) - k_dt_date_to_count(today))/365.0 ));

	vector<date> payment;
	int *pdate;
	int num_rol;
	int i;
	int num_month;
	double fwdrate;
	

	int ssdate;
	date sw_sdate;

    int spdate;
	date spot_date;

	if (ptype == M)
		num_month = 1;
	else if (ptype == Q)
		num_month = 3;
	else if (ptype == S)
		num_month = 6;
	else if (ptype == A)
		num_month = 12;
	else 
		return(-2000.0); /// payment type error

	date *fixing;
	int  *fdate;

	sw_sdate = k_dt_getSpotDate(omat, pay_cal);   // swap start
	ssdate = k_dt_date_to_count(sw_sdate) - k_dt_date_to_count(today);

	spot_date = k_dt_getSpotDate(today,pay_cal);  // spot day
	spdate = k_dt_date_to_count(spot_date) - k_dt_date_to_count(today);

	

	num_rol = ( (smat.year-sw_sdate.year)*12 + (smat.month-sw_sdate.month))/num_month;
    payment = k_swap_scheduleGenerator(sw_sdate,num_rol,num_month,BizDayConv,pay_cal);
	pdate = (int *)malloc(sizeof(int)*num_rol);
	
	for(i=0; i < num_rol; i++){
		pdate[i] = k_dt_date_to_count(payment[i]) - k_dt_date_to_count(today);	
	}

	
	double annuity;

	annuity = (*zCurve).DF(ssdate/365.0, pdate[0]/365.0)*k_dt_yearFrac(sw_sdate,payment[0],DayCountConv);

//	printf("%.15f %.15f \n",(*zCurve).DF(ssdate/365.0, pdate[0]/365.0),k_dt_yearFrac(sw_sdate,payment[0],DayCountConv));
	for(i=1; i < num_rol; i++){

		annuity +=(*zCurve).DF(ssdate/365.0, pdate[i]/365.0)*k_dt_yearFrac(payment[i-1],payment[i],DayCountConv);
//		printf("%.15f %.15f\n",(*zCurve).DF(ssdate/365.0, pdate[i]/365.0),k_dt_yearFrac(payment[i-1],payment[i],DayCountConv) );
	}



	date fstart;
	vector<date> fraed;

	int fsdate, frdate;
    double fra;
	double FLEG;
    


	if (frtype == M3)
		num_month = 3;
	else if (frtype == M6)
		num_month = 6;
	else if (frtype == M1)
		num_month = 1;

	int num_rol1;
    vector<date> payment1;
	int *pdate1;

	num_rol1 = ( (smat.year-sw_sdate.year)*12 + (smat.month-sw_sdate.month))/num_month;
	payment1 = k_swap_scheduleGenerator(sw_sdate,num_rol1,num_month,BizDayConv,pay_cal);
	pdate1 = (int *)malloc(sizeof(int)*num_rol1);
	fixing = (date *)malloc(sizeof(date)*num_rol1);
	fdate = (int *)malloc(sizeof(int)*num_rol1);

	for(i=0; i < num_rol1; i++){
		pdate1[i] = k_dt_date_to_count(payment1[i]) - k_dt_date_to_count(today);	
	}

	fixing[0] = omat;
	fdate[0] = k_dt_date_to_count(fixing[0]) - k_dt_date_to_count(today);
	for(i=1; i < num_rol1; i++){
		fixing[i] = k_swap_fix_date_shift(payment1[i-1],fix_cal);
		fdate[i] = k_dt_date_to_count(fixing[i]) - k_dt_date_to_count(today);
	}


	fstart = k_dt_getSpotDate(fixing[0],pay_cal);
	fraed = k_dt_scheduleGenerator(fstart,1,num_month,BizDayConv,pay_cal);

	fsdate = k_dt_date_to_count(fstart) - k_dt_date_to_count(today);
	frdate = k_dt_date_to_count(fraed[0]) - k_dt_date_to_count(today);

	fra = (*zCurve).OR(fsdate/365.0,frdate/365.0, k_dt_yearFrac(fstart,fraed[0],DayCountConv));

	FLEG = fra* k_dt_yearFrac(sw_sdate,payment1[0],DayCountConv)*(*zCurve).DF(ssdate/365.0, pdate1[0]/365.0);

//	printf("%.15f %d %d \n", fra*100.0,k_dt_date_to_count(fixing[0]), k_dt_date_to_count(payment1[0]),(*zCurve).DF(ssdate/365.0, pdate1[0]/365.0));

	

	for (i= 1; i< num_rol1; i++) {

		fstart = k_dt_getSpotDate(fixing[i],pay_cal);
		fraed = k_dt_scheduleGenerator(fstart,1,num_month,BizDayConv,pay_cal);

		fsdate = k_dt_date_to_count(fstart) - k_dt_date_to_count(today);
		frdate = k_dt_date_to_count(fraed[0]) - k_dt_date_to_count(today);

        fra = (*zCurve).OR(fsdate/365.0,frdate/365.0, k_dt_yearFrac(fstart,fraed[0],DayCountConv));

        FLEG += fra * k_dt_yearFrac(payment1[i-1],payment1[i],DayCountConv)*(*zCurve).DF(ssdate/365.0, pdate1[i]/365.0);

	//	printf("%.15f %d %d %.15f\n", fra*100.0,k_dt_date_to_count(fixing[i]), k_dt_date_to_count(payment1[i]),(*zCurve).DF(ssdate/365.0, pdate1[i]/365.0));

	}


	fwdrate = FLEG/annuity;

	
	if (pdate) free(pdate);
	if (fixing) free(fixing);
	if (fdate) free(fdate);
	if (pdate1) free(pdate1);

	return(fwdrate);
}



double k_ForwardSwapRate( int n, double *Tau, double *ZB)
{
	//	generate forward swap rate and spot swap rate from array of zero coupon bond price
	//	formular derived by hangseob	:	need to check!!!
	//  
	//	n:	number of payment
	//	Tau[0] ~ Tau[n-1]: payment interval size considering day count convention
	//                     Tau[0] = first payment day - swap starting day
	//                     Tau[1] = second payment day - first payment day
	//                 ... Tau[n-1] = n-th payment day - (n-1)-th payment day 
	//	ZB[k]: zero coupon bond price maturing at the (k + 1)-th payment day

	double annuity = 0.0;
	for (int k= 0;k< n;k++) 
		annuity  += ZB[k]*Tau[k];
	return (ZB[0] - ZB[n])/annuity;
}


double k_ForwardSwapRate( int n, double tau, double *ZB)
{
	//	generate forward swap rate and spot swap rate from array of zero coupon bond price
	//	formular derived by hangseob	:	need to check!!!
	//  
	//	n:	number of payment
	//	tau: fixed payment interval, namely 0.25 = 3M, 0.5 = 6M, ....
	//	ZB[k]: zero coupon bond price maturing at the (k + 1)-th payment day

	double annuity = 0.0;
	for (int k= 0;k< n;k++) 
		annuity  += ZB[k]*tau;
	return (ZB[0] - ZB[n])/annuity;
}

double k_ForwardSwapRate(k_cZeroCurve *zCurve, int n, double Ts, PayType ptype)
{
	// simplified version;
	//	generate forward swap rate and spot swap rate from zero curve
	//	formular derived by hangseob	:	need to check!!!	modified by Jaemin Ahn
	//	n:	number of payment
	//	Ts : swap starting time(i.e. first fixing point)
	//		*caution!!: when it is fwd swap, Ts is not zero!!!
	//  ptype : Payment Type(fixed coupon) M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually


	double tau;

	if (ptype == M)
		tau = 1.0/12;
	else if (ptype == Q)
		tau = 0.25;
	else if (ptype == S)
		tau = 0.5;
	else if (ptype == A)
		tau = 1.0;
	else 
		return(-2000.0); /// payment type error


	double annuity = 0.0;
	for (int k= 1;k<= n;k++) 
		annuity  += (*zCurve).DF(Ts+k*tau)*tau;
	return ((*zCurve).DF(Ts) - (*zCurve).DF(Ts+tau*n))/annuity;
}

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
				   )
{

	//	if option holder pays fixed, Flag = 0;
	//	if option holder receives fixed, Flag = 1;
	//	strike :	strike of swaption
	//  sigma : volatility 
	//  L : Notional amount
	//  *zCurve : zero curve structure
	//	today :	today
	//  sdate : swap starting date (effective date for swap)
	//  smat  : swap maturity 
	//  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually
	//  BizDayConv : 0=modified following,  1=following,  2=preceding
	//	DayCountConv :  0=A/365, 1=30/360, 2=A/360, 3=A/A
	//  pay_cal : calendar --> SEOUL, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, SEOUL_NEW_YORK,
	//                         TARGET_NEW_YORK, TOKYO_NEW_YORK

	double s0;
	
	s0 = k_ForwardSwapRate(zCurve,today,omat,smat,ptype,frtype,BizDayConv,DayCountConv,pay_cal,fix_cal);


	vector<date> payment;
	int *pdate;
	int num_rol;
	int i;
	int num_month;
	

	int ssdate;
	date sw_sdate;

	int spdate;
	date spot_date;

	if (ptype == M)
		num_month = 1;
	else if (ptype == Q)
		num_month = 3;
	else if (ptype == S)
		num_month = 6;
	else if (ptype == A)
		num_month = 12;
	else 
		return(-2000.0); /// payment type error

	
	sw_sdate = k_dt_getSpotDate(omat, pay_cal);   // swap start
	ssdate = k_dt_date_to_count(sw_sdate) - k_dt_date_to_count(today);

	spot_date = k_dt_getSpotDate(today,pay_cal);  // spot day
	spdate = k_dt_date_to_count(spot_date) - k_dt_date_to_count(today);



	num_rol = ( (smat.year-sw_sdate.year)*12 + (smat.month-sw_sdate.month))/num_month;
	payment = k_swap_scheduleGenerator(sw_sdate,num_rol,num_month,BizDayConv,pay_cal);
	pdate = (int *)malloc(sizeof(int)*num_rol);

	for(i=0; i < num_rol; i++){
		pdate[i] = k_dt_date_to_count(payment[i]) - k_dt_date_to_count(today);	
	}


	double annuity;

	annuity = (*zCurve).DF(spdate/365.0, pdate[0]/365.0)*k_dt_yearFrac(sw_sdate,payment[0],DayCountConv);

	//	printf("%.15f %.15f \n",(*zCurve).DF(ssdate/365.0, pdate[0]/365.0),k_dt_yearFrac(sw_sdate,payment[0],DayCountConv));
	for(i=1; i < num_rol; i++){

		annuity +=(*zCurve).DF(spdate/365.0, pdate[i]/365.0)*k_dt_yearFrac(payment[i-1],payment[i],DayCountConv);
		//		printf("%.15f %.15f\n",(*zCurve).DF(ssdate/365.0, pdate[i]/365.0),k_dt_yearFrac(payment[i-1],payment[i],DayCountConv) );
	}

	if (pdate) free(pdate);

	if (Flag == 0)	
		// payer swaption
		return L*k_BlackCall( annuity, s0, strike, (k_dt_date_to_count(omat)-k_dt_date_to_count(today))/365.0, sigma);
	else
		// receiver swaption
		return L*k_BlackPut( annuity, s0, strike, (k_dt_date_to_count(omat)-k_dt_date_to_count(today))/365.0, sigma);
}


double k_Swaption( int Flag, k_cZeroCurve *zCurve, int n, double Ts, PayType	ptype, double strike, double sigma, double L)  
 {
   //   simplified swaption pricer
   //	if option holder pays fixed, Flag = 0;
   //	if option holder receives fixed, Flag = 1;
   //	n:	number of payment
   //	Ts :	swap starting time(i.e option maturity time)
   //   ptype: fixed rae payement type M=1M, Q=3M, S=6M, A=12M
   //   strike: swaption strike
   //   sigma : volatility
   //   L : principal
   	 

	 double tau;

	 if (ptype == M)
		 tau = 1.0/12;
	 else if (ptype == Q)
		 tau = 0.25;
	 else if (ptype == S)
		 tau = 0.5;
	 else if (ptype == A)
		 tau = 1.0;
	 else 
		 return(-2000.0); /// payment type error

   double s0 = k_ForwardSwapRate(zCurve, n, Ts, ptype);

   double annuity = 0.0;

   for (int k= 1;k<= n;k++) 
	   annuity  += (*zCurve).DF(Ts+k*tau)*tau;	

	if (Flag == 0)	
	   return L*k_BlackCall( annuity, s0, strike, Ts, sigma);
	else
	   return L*k_BlackPut( annuity, s0, strike, Ts, sigma);
 }


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
						  )

	//	*zCurve : zero curve structure, 	
	//	strike : cap strike 
	//	today :	today (trade date)
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

{
	vector<date> payment;
	int *pdate;
	int *fdate;
	date *fixing;
	int num_rol;
	int i;
	int num_month;
	double *caplet_fwd;
	double *caplet_yfrac;
	double *caplet_dfac;



	if (ptype == M)
		num_month = 1;
	else if (ptype == Q)
		num_month = 3;
	else if (ptype == S)
		num_month = 6;
	else if (ptype == A)
		num_month = 12;
	else 
		return(-2000.0); /// payment type error

	num_rol = ( (mk_cap_mat[num_cap-1].year-today.year)*12 + (mk_cap_mat[num_cap-1].month-today.month))/num_month ;

	date effective_date;

	int edate;
	int *num_tenor;

	

	effective_date = k_swap_getSpotDate(today,pay_cal);
	edate = k_dt_date_to_count(effective_date) - k_dt_date_to_count(today);
	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,pay_cal);

	pdate = (int *)malloc(sizeof(int)*num_rol);   // payment date = calculation end 
	fdate = (int *)malloc(sizeof(int)*num_rol);   // fixing date 
	fixing = (date *)malloc(sizeof(date)*num_rol);
	caplet_dfac = (double *)malloc(sizeof(double )*num_caplet);
	caplet_fwd = (double *)malloc(sizeof(double )*num_caplet);
	caplet_yfrac = (double *)malloc(sizeof(double )*num_caplet);
	

    double *mk_cap_val;

	mk_cap_val = (double *)malloc(sizeof(double) * num_cap);
    num_tenor = (int *)malloc(sizeof(int )*num_cap);

	for(i=0; i < num_cap; i++){
		mk_cap_val[i] = k_cap(zCurve,strike, 1.0, mk_cap_vol[i], today, mk_cap_mat[i], ptype, BizDayConv,
		                      DayCountConv, pay_cal, fix_cal);
		num_tenor[i] =  ( (mk_cap_mat[i].year-today.year)*12 + (mk_cap_mat[i].month-today.month))/num_month 
 - 1;

	//	printf("cap %.5f \n",mk_cap_val[i]*100);
	}


	for(i=0; i < num_rol; i++){
		pdate[i] = k_dt_date_to_count(payment[i]) - k_dt_date_to_count(today);
		fixing[i] = k_swap_fix_date_shift(payment[i],fix_cal);
		fdate[i] = k_dt_date_to_count(fixing[i]) - k_dt_date_to_count(today);
	}

	for(i=0; i < num_caplet; i++){
        caplet_mat[i] = payment[i+1];
	}

	for(i=0; i < num_tenor[0]; i++){
       caplet_vol[i] = mk_cap_vol[0];
	}


	int frdate;

	vector<date> fraed;
	date fstart;
	int fsdate;

	for (i= 1; i< num_rol; i++) {
		fstart = k_dt_getSpotDate(fixing[i-1],pay_cal);
		fraed = k_dt_scheduleGenerator(fstart,1,num_month,BizDayConv,pay_cal);

		fsdate = k_dt_date_to_count(fstart) - k_dt_date_to_count(today);
		frdate = k_dt_date_to_count(fraed[0]) - k_dt_date_to_count(today);

		caplet_dfac[i-1] = (*zCurve).DF(edate/365.0,pdate[i]/365.0);
		caplet_fwd[i-1] = (*zCurve).OR(fsdate/365.0,frdate/365.0, k_dt_yearFrac(fstart,fraed[0],DayCountConv));
        caplet_yfrac[i-1] = k_dt_yearFrac(payment[i-1],payment[i],DayCountConv);
	}




	int k, n, j;


    int temp_capn;
	int idx;
	    
	double dis0, dis1 , dis;
	double eps = 1.0e-05;
	double TOL = 1.0e-12;
	double p0, p1, q0, q1, p, q;
	
	int temp_days;
	double start;
	
	double diff;
	int maxI = 300;
	double TF = 1.0;

	for(j=0; j < num_cap - 1; j++){
         temp_capn = num_tenor[j+1] - num_tenor[j];
		 temp_days = k_dt_numDays(caplet_mat[num_tenor[j]-1], caplet_mat[num_tenor[j+1]-1], DayCountConv);

         start = caplet_vol[num_tenor[j]-1];

		 p1 = 2.0;
		 p0 = 0.00000001;

		 diff = mk_cap_val[j+1] - mk_cap_val[j];
		
		 q0 = q1 = 0.0;

		 dis0 = dis1 = 0.0;
		 for(k=0; k < temp_capn; k++){

			
			 idx = num_tenor[j] + k;

			 dis0 += (p0-start)*k_dt_numDays(caplet_mat[idx-1],caplet_mat[idx],DayCountConv)*1.0/temp_days;
			 dis1 += (p1-start)*k_dt_numDays(caplet_mat[idx-1],caplet_mat[idx],DayCountConv)*1.0/temp_days;

			// printf("id %d %.5f %.5f\n",idx, start+dis0, start+dis1);
			 q0 += k_caplet(caplet_dfac[idx],caplet_fwd[idx],strike,fdate[idx]/365.0,caplet_yfrac[idx],1.0,
				            start + dis0);
			 q1 += k_caplet(caplet_dfac[idx],caplet_fwd[idx],strike,fdate[idx]/365.0,caplet_yfrac[idx],1.0,
				 start + dis1);
		 }
		 q0 -= diff; q1 -= diff;

		//  printf("q %.5f %.5f \n",q0,q1);

		 if (q0*q1 >= 0.0){

			 for(k=0; k < temp_capn; k++){
				 idx = num_tenor[j] + k;
				 caplet_vol[idx] = start;
			 }
			 TF = 0.0;
		 }
		 else{

			for( n=1; n < maxI; n++){
                  
				 p  = p1 - q1*(p1-p0)/(q1-q0);
				 q = 0.0;

				 dis = 0.0;
				 
				 for(k=0; k < temp_capn; k++){

					 idx = num_tenor[j] + k;

					 dis += (p-start)*k_dt_numDays(caplet_mat[idx-1],caplet_mat[idx],DayCountConv)*1.0/temp_days;
					 
					 q += k_caplet(caplet_dfac[idx],caplet_fwd[idx],strike,fdate[idx]/365.0,caplet_yfrac[idx],1.0,
						 start + dis);
					
				 }

				 q -= diff;

				 if ( fabs(q) < TOL)  break;

				 if (q*q1 < 0.0){
					  p0 = p;
					  q0 = q;
				 }
				 else{
					  q1 = q;
					  p1 = p;
				 }

			 } // end of n

			 dis1 = 0.0;

			 for(k=0; k < temp_capn; k++){

				 idx = num_tenor[j] + k;

				  

				 dis1 += (p-start)*k_dt_numDays(caplet_mat[idx-1],caplet_mat[idx],DayCountConv)*1.0/temp_days;

				 caplet_vol[idx] = start + dis1;				

			 }	 


		 }
	}


	for(i=0; i < num_caplet; i++){
		caplet_mat[i] = payment[i];
	}

	if (pdate) free(pdate);
	if (fdate) free(fdate);
	if (fixing) free(fixing);
	if (caplet_fwd) free(caplet_fwd);
	if (caplet_dfac) free(caplet_dfac);
	if (caplet_yfrac) free(caplet_yfrac);
	if (mk_cap_val) free(mk_cap_val);
    if (num_tenor) free(num_tenor);

	return(TF);

}


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
						  )

	  // simplified caplet_stripping 
	  //	the caplet whose payoff fixed at today is not included
	  //	*zCurve : zero curve structure, 	
	  //	strike : cap strike 
	  //  *mk_cap_vol : market cap volatilities w.r.t. maturities 
	  //  *mk_cap_mat : cap maturities
	  //  num_cap : number of caps (number of maturities) 
	  //  ptype : Payment Type M: Monthly, Q: Quarterly, S: Semi-Annually, A: Annually
	  //  num_caplet : number of caplets 
	  //  Output:*caplet_mat = caplet maturities (from today(time 0) to caplet maturity)
	  //          *caplet_vol = stripped caplet volatilities 

{
	double tau;

	if (ptype == M)
		tau = 1.0/12;
	else if (ptype == Q)
		tau = 0.25;
	else if (ptype == S)
		tau = 0.5;
	else if (ptype == A)
		tau = 1.0;


	int num_rol;
	int i;
	double *caplet_fwd;
	double *caplet_dfac;

	
	num_rol = num_caplet + 1;

	
	int *num_tenor;

	caplet_dfac = (double *)malloc(sizeof(double )*num_caplet);
	caplet_fwd = (double *)malloc(sizeof(double )*num_caplet);


	double *mk_cap_val;

	mk_cap_val = (double *)malloc(sizeof(double) * num_cap);
	num_tenor = (int *)malloc(sizeof(int )*num_cap);

	for(i=0; i < num_cap; i++){

		num_tenor[i] = (int )(mk_cap_mat[i]/tau + 0.001) - 1;

		mk_cap_val[i] = k_cap(zCurve,strike,num_tenor[i],1.0,mk_cap_vol[i],ptype);

	//	printf("%d %d %.5f \n",i,num_tenor[i],mk_cap_val[i]);
		
	}

	for(i=0; i < num_caplet; i++){
		caplet_mat[i] = tau*(i+1.0);
	}

	for(i=0; i < num_tenor[0]; i++){
		caplet_vol[i] = mk_cap_vol[0];
	}


	for (i= 1; i< num_rol; i++) {
		caplet_dfac[i-1] = (*zCurve).DF(tau*(i+1.0));
		caplet_fwd[i-1] = (*zCurve).OR(tau*i,tau*(i+1.0),tau);

		// printf("%d %.5f %.5f \n",i-1,caplet_dfac[i-1],caplet_fwd[i-1]);
	}

	int k, n, j;


	int temp_capn;
	int idx;
	

	double dis0, dis1 , dis;
	double eps = 1.0e-05;
	double TOL = 1.0e-12;
	double p0, p1, q0, q1, p, q;

	double start;

	double diff;
	int maxI = 300;
	double TF = 1.0;

	for(j=0; j < num_cap - 1; j++){
		temp_capn = num_tenor[j+1] - num_tenor[j];
		
		start = caplet_vol[num_tenor[j]-1];

		p1 = 2.0;
		p0 = 0.00000001;

		diff = mk_cap_val[j+1] - mk_cap_val[j];

		q0 = q1 = 0.0;

		dis0 = (p0-start)/temp_capn;
		dis1 = (p1-start)/temp_capn;

		for(k=0; k < temp_capn; k++){


			idx = num_tenor[j] + k;

			q0 += k_caplet(caplet_dfac[idx],caplet_fwd[idx],strike,caplet_mat[idx],tau,1.0,start+(k+1.0)*dis0);
            q1 += k_caplet(caplet_dfac[idx],caplet_fwd[idx],strike,caplet_mat[idx],tau,1.0,start+(k+1.0)*dis1);
		//	printf("id %d %.5f  %.5f %.5f\n",idx,start*100.0, start+(k+1.0)*dis0, start+(k+1.0)*dis1);
		 }
		 q0 -= diff; q1 -= diff;

		// printf("q %.5f %.5f \n",q0,q1);

		if (q0*q1 >= 0.0){

			for(k=0; k < temp_capn; k++){
				idx = num_tenor[j] + k;
				caplet_vol[idx] = start;
			}
			TF = 0.0;
		}
		else{

			for( n=1; n < maxI; n++){

				p  = p1 - q1*(p1-p0)/(q1-q0);
			
				q = 0.0;

				dis = (p-start)/temp_capn;

				for(k=0; k < temp_capn; k++){

					idx = num_tenor[j] + k;

					q += k_caplet(caplet_dfac[idx],caplet_fwd[idx],strike,caplet_mat[idx],tau,1.0,start+(k+1.0)*dis);

				}

				q -= diff;

				if ( fabs(q) < TOL)  break;

				if (q*q1 < 0.0){
					p0 = p;
					q0 = q;
				}
				else{
					q1 = q;
					p1 = p;
				}

			} // end of n

			dis = (p-start)/temp_capn;

			for(k=0; k < temp_capn; k++){

				idx = num_tenor[j] + k;

				caplet_vol[idx] = start + (k+1.0)*dis;				

			}	 


		}
	}

	if (caplet_fwd) free(caplet_fwd);
	if (caplet_dfac) free(caplet_dfac);
	if (mk_cap_val) free(mk_cap_val);
	if (num_tenor) free(num_tenor);

	return(TF);
}