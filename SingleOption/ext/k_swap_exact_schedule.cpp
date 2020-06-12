#include <kdb_swap_exact_schedule.hpp>

int k_swap_holiday_search(vector<long> hdays,int size,int key)
{     int lb=0,ub=size-1,mid;             //lb=>lower bound,ub=>upper bound
int tf=0; 

for(;lb<= ub;)
{
	mid=(lb+ub)/2;

	if(hdays[mid] == key)
	{
		tf = 1;  
		return(tf);

	} 
	else
		if(hdays[mid]< key)
			lb=mid+1;
		else
			if (hdays[mid]>key)
				ub=mid-1;
}

return(tf);

}


long k_swap_day_shift(date dt,calendar cal,int n, int direction){
	// dt: business date !!!
	// direction=-1: proceeding, direction = 1: following

	long count, adj_count;

	int ilast;
	int j;
	int NH;

	count=k_dt_date_to_count(dt);
	adj_count=count;
	vector<long> hdays;

	k_dt_getHolidays(hdays,cal);
	ilast= (int)hdays.size()-1;

	for(j=0; j < n; j++){

		adj_count = adj_count + direction;

		NH = k_swap_holiday_search(hdays,ilast,adj_count);

		for(;(NH != 0) || (k_dt_week(adj_count) == 0) || (k_dt_week(adj_count) == 1);){
			if (NH != 0){
				adj_count = adj_count + direction;
			}
			if (k_dt_week(adj_count) == 0){
				if (direction == 1){
					adj_count = adj_count + 2*direction;
				}
				else{
					adj_count = adj_count + direction;
				}
			}
			if (k_dt_week(adj_count) == 1){
				if (direction == 1){
					adj_count = adj_count + direction;
				}
				else{
					adj_count = adj_count+ 2*direction;
				}
			}	
			NH = k_swap_holiday_search(hdays,ilast,adj_count);
		}
	}
	return(adj_count); 
}


/************************************************************************/
// shifting date to business date according to the business day convention
// BizDayConv: 1=following, 2=preceding, default=modified following
/************************************************************************/

date k_swap_dt_workday(date dt, int BizDayConv, calendar calendar_)
{
	long adj_count;
	long count;
	long NH;
	long direction;
	int ilast;

	vector<long> hdays;

	k_dt_getHolidays(hdays,calendar_);
	ilast= (int)hdays.size()-1;

	count = k_dt_date_to_count(dt);
	adj_count = count;



	switch(BizDayConv)
	{
	case 1: // following
		direction = 1;
	case 2: //preceding
		direction = -1;
	default: // modified following
		direction = 1;	
	}

	NH = k_swap_holiday_search(hdays,ilast,adj_count);
	for(;(NH != 0) || (k_dt_week(adj_count) == 0) || (k_dt_week(adj_count) == 1);){
		if (NH != 0){
			adj_count = adj_count + direction;
		}
		if (k_dt_week(adj_count) == 0){
			if (direction == 1){
				adj_count = adj_count + 2*direction;
			}
			else{
				adj_count = adj_count + direction;
			}
		}
		if (k_dt_week(adj_count) == 1){
			if (direction == 1){
				adj_count = adj_count + direction;
			}
			else{
				adj_count = adj_count+ 2*direction;
			}
		}	
		NH = k_swap_holiday_search(hdays,ilast,adj_count);
	}

	date adj_dt;
	adj_dt=k_dt_count_to_date(adj_count);

	if ((BizDayConv != 1) && (BizDayConv != 2)){ 

		if(dt.month==adj_dt.month)
			return adj_dt;
		else return( k_swap_dt_workday(dt, 2, calendar_) );
	}
	else{
		return(adj_dt);
	}

}


vector<date> k_swap_scheduleGenerator(date startDate,int numRollover,int schedule_term,int BizDayConv,calendar calendar_)
{
	int i;
	vector<date> dt(numRollover);
	vector<date> adj_dt(numRollover);
	for(i=0;i<numRollover;i++)
	{
		dt[i]=k_dt_month_adder(startDate,schedule_term*(i+1));
		adj_dt[i]=k_swap_dt_workday(dt[i],BizDayConv,calendar_);
	}
	return adj_dt;
}


date k_swap_getSpotDate(date dt, calendar calendar_)
{
	date r_dt;
	switch(calendar_)
	{
	case SEOUL:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,1,1));
		break;
	case NEW_YORK:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,1));
		break;
	case TARGET:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,1));
		break;
	case TOKYO:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,1));
		break;
	case NEW_YORK_LONDON:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,1));
		break;
	case LONDON:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,1));
		break;
	case SEOUL_NEW_YORK: // for CRS
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,1));
		break;
	}

	return r_dt;
}


date k_swap_fix_date_shift(date dt, calendar calendar_)
// determine fixing day shifted back from payment date 
{
	date r_dt;
	switch(calendar_)
	{
	case SEOUL:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,1,-1));
		break;
	case NEW_YORK:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,-1));
		break;
	case TARGET:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,-1));
		break;
	case TOKYO:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,-1));
		break;
	case NEW_YORK_LONDON:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,-1));
		break;
	case LONDON:
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,-1));
		break;
	case SEOUL_NEW_YORK: // for CRS
		r_dt=k_dt_count_to_date(k_swap_day_shift(dt,calendar_,2,-1));
		break;
	}

	return r_dt;
}


double k_swap_KRWIRS_FLOATING_LEG(
								  date pricing_date,     // 가격산정일
								  date trade_date,       // trade date
								  date maturity,          // swap 총 만기  unit: year
								  FRateType  frtype,     // Floating Rate Type (M3=3 month, M6=6 month)
								  double pfixed_rate,    // rate fixed previously
								  double margin,         // margin
								  int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
								  int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
								  k_cZeroCurve *ZC  // zero curve strucutre at pricing date
								  )
{
	vector<date> payment;  // vector for payment schedule
	int *ifd, *ipd;        // initial fixing schedule, payment schedule(거래 입력당시 확정된 스케쥴)
	int *cpd;              // current payment day-interval from pricing date (현재 남은 payment dates)
	int *csd, *ced;        // calculation start dates, end dates(거래 입력당시 확정되는 calculation period)
	int num_rol;           // initial number of rollovers(거래 입력 당시 확정되는 총 rollover 수)
	int cur_rol;           // current number of remaining rollovers(현 가격산정 시점에 남아있는 rollover 수)
	int cur_id;            // current index 
	int i;
	double sum;
	date spot_day;  // spot day 
	int spd;        // spot day interval from pricing day	
	int num_month;
	date effective_date;

	if (frtype == M3)
		num_month = 3;
	else if (frtype == M6)
		num_month = 6;
	else if (frtype == M1)
		num_month = 1;



	num_rol = ((maturity.year-trade_date.year)*12+(maturity.month-trade_date.month))/num_month;
	effective_date = k_swap_getSpotDate(trade_date,SEOUL);
	spot_day = k_swap_getSpotDate(pricing_date,SEOUL);
	spd = k_dt_date_to_count(spot_day)-k_dt_date_to_count(pricing_date);

	ifd = (int *)malloc(sizeof(int)*num_rol);
	ipd = (int *)malloc(sizeof(int)*num_rol);
	csd = (int *)malloc(sizeof(int)*num_rol);
	ced = (int *)malloc(sizeof(int)*num_rol);

	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,SEOUL);

	ipd[0] = k_dt_date_to_count(payment[0]);
	ifd[0] = k_dt_date_to_count(trade_date);
	csd[0] = k_dt_date_to_count(effective_date);
	ced[0] = ipd[0];
	i=0;

	for(i=1; i < num_rol; i++){
		ipd[i] = k_dt_date_to_count(payment[i]);
		ifd[i] = k_swap_day_shift(payment[i-1],SEOUL,1,-1);
		csd[i] = ipd[i-1];
		ced[i] = ipd[i];
	}

	for(i=0; i < num_rol; i++){

		if (k_dt_date_to_count(payment[i]) >= k_dt_date_to_count(pricing_date)) {
			cur_rol = num_rol-i;
			cur_id = i;
			break;
		}
	}


	cpd = (int *)malloc(sizeof(int)*cur_rol);

	for(i=0; i < cur_rol; i++){
		cpd[i] = ipd[i+cur_id] - k_dt_date_to_count(pricing_date);
	}


	double fra;
	vector<date> fraed;

	sum = (pfixed_rate+margin) * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv)* (*ZC).DF(spd/365.0,cpd[0]/365.);


	for(i=1; i < cur_rol; i++){

		fraed = k_swap_scheduleGenerator(k_dt_getSpotDate(k_dt_count_to_date(ifd[i+cur_id]),SEOUL),1,num_month,BizDayConv,SEOUL);

		fra = (*ZC).OR(  k_dt_yearFrac(pricing_date,k_dt_getSpotDate(k_dt_count_to_date(ifd[i+cur_id]),SEOUL),  DayCountConv),
			k_dt_yearFrac(pricing_date,fraed[0],  DayCountConv),
			k_dt_yearFrac(k_dt_getSpotDate(k_dt_count_to_date(ifd[i+cur_id]),SEOUL),fraed[0],DayCountConv));

		sum += (fra+margin)*k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)
			*(*ZC).DF(spd/365.0,cpd[i]/365.);
	}


	if (ifd) free(ifd);
	if (ipd) free(ipd);
	if (cpd) free(cpd);
	if (csd) free(csd);
	if (ced) free(ced);

	return(sum);

}

double k_swap_KRWIRS_FIXED_LEG(
							   date pricing_date,     // 가격산정일
							   date trade_date,       // trade date
							   date  maturity,         // swap 총 만기 unit:year
							   PayType ptype,         // 이자지급주기 M=montly, Q=quarterly, S=semi-annually, A=annually 
							   double pfixed_rate,    // 확정고정금리
							   int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
							   int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
							   k_cZeroCurve *ZC  // zero curve strucutre at pricing date
							   )
{
	vector<date> payment;  // vector for payment schedule
	int *ipd;              // initial payment schedule(거래 입력당시 확정된 스케쥴)
	int *cpd;              // current payment day-interval from pricing date (현재 남은 payment dates)
	int *csd, *ced;        // calculation start dates, end dates(거래 입력당시 확정되는 calculation period)
	int num_rol;           // initial number of rollovers(거래 입력 당시 확정되는 총 rollover 수)
	int cur_rol;           // current number of remaining rollovers(현 가격산정 시점에 남아있는 rollover 수)
	int cur_id;            // current index 
	int i;
	double sum;
	date spot_day;  // spot day 
	int spd;        // spot day interval from pricing day	
	int num_month;
	date effective_date;

	if (ptype == M)
		num_month = 1;
	else if (ptype == Q)
		num_month = 3;
	else if (ptype == S)
		num_month = 6;
	else if (ptype == A)
		num_month = 12;

	// num_rol = (12/num_month)*maturity;

	num_rol = ((maturity.year-trade_date.year)*12+(maturity.month-trade_date.month))/num_month;

	effective_date = k_swap_getSpotDate(trade_date,SEOUL);
	spot_day = k_swap_getSpotDate(pricing_date,SEOUL);
	spd = k_dt_date_to_count(spot_day)-k_dt_date_to_count(pricing_date);

	ipd = (int *)malloc(sizeof(int)*num_rol);
	csd = (int *)malloc(sizeof(int)*num_rol);
	ced = (int *)malloc(sizeof(int)*num_rol);

	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,SEOUL);

	ipd[0] = k_dt_date_to_count(payment[0]);
	csd[0] = k_dt_date_to_count(effective_date);
	ced[0] = ipd[0];
	for(i=1; i < num_rol; i++){
		ipd[i] = k_dt_date_to_count(payment[i]);
		csd[i] = ipd[i-1];
		ced[i] = ipd[i];
	}


	for(i=0; i < num_rol; i++){

		if (k_dt_date_to_count(payment[i]) >= k_dt_date_to_count(pricing_date)) {
			cur_rol = num_rol-i;
			cur_id = i;
			break;
		}
	}

	cpd = (int *)malloc(sizeof(int)*cur_rol);

	for(i=0; i < cur_rol; i++){
		cpd[i] = ipd[i+cur_id] - k_dt_date_to_count(pricing_date);
	}

	sum = pfixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv)* (*ZC).DF(spd/365.0,cpd[0]/365.);
	for(i=1; i < cur_rol; i++){
		sum += pfixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)*(*ZC).DF(spd/365.0,cpd[i]/365.);
	}


	if (ipd) free(ipd);
	if (cpd) free(cpd);
	if (csd) free(csd);
	if (ced) free(ced);

	return(sum);

}

double k_swap_USDIRS_FIXED_LEG(
							   date pricing_date,     // 가격산정일
							   date trade_date,       // trade date
							   date  maturity,         // swap 총 만기 unit:year
							   PayType ptype,         // 이자지급주기 M=montly, Q=quarterly, S=semi-annually, A=annually 
							   double pfixed_rate,    // 확정고정금리
							   int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
							   int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
							   k_cZeroCurve	*ZC // zero curve strucutre at pricing date
							   )
{
	vector<date> payment;  // vector for payment schedule
	int *ipd;              // initial payment schedule(거래 입력당시 확정된 스케쥴)
	int *cpd;              // current payment day-interval from pricing date (현재 남은 payment dates)
	int *csd, *ced;        // calculation start dates, end dates(거래 입력당시 확정되는 calculation period)
	int num_rol;           // initial number of rollovers(거래 입력 당시 확정되는 총 rollover 수)
	int cur_rol;           // current number of remaining rollovers(현 가격산정 시점에 남아있는 rollover 수)
	int cur_id;            // current index 
	int i;
	double sum;
	date spot_day;  // spot day 
	int spd;        // spot day interval from pricing day	
	int num_month;
	date effective_date;   // effective date

	if (ptype == M)
		num_month = 1;
	else if (ptype == Q)
		num_month = 3;
	else if (ptype == S)
		num_month = 6;
	else if (ptype == A)
		num_month = 12;


	num_rol = ((maturity.year-trade_date.year)*12+(maturity.month-trade_date.month))/num_month;

	effective_date = k_swap_getSpotDate(trade_date,NEW_YORK_LONDON);
	spot_day = k_swap_getSpotDate(pricing_date,NEW_YORK_LONDON);
	spd = k_dt_date_to_count(spot_day)-k_dt_date_to_count(pricing_date);

	ipd = (int *)malloc(sizeof(int)*num_rol);
	csd = (int *)malloc(sizeof(int)*num_rol);
	ced = (int *)malloc(sizeof(int)*num_rol);

	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,NEW_YORK_LONDON);

	ipd[0] = k_dt_date_to_count(payment[0]);
	csd[0] = k_dt_date_to_count(effective_date);
	ced[0] = ipd[0];
	for(i=1; i < num_rol; i++){
		ipd[i] = k_dt_date_to_count(payment[i]);
		csd[i] = ipd[i-1];
		ced[i] = ipd[i];
	}

	for(i=0; i < num_rol; i++){

		if (k_dt_date_to_count(payment[i]) >= k_dt_date_to_count(pricing_date)) {
			cur_rol = num_rol-i;
			cur_id = i;
			break;
		}
	}

	cpd = (int *)malloc(sizeof(int)*cur_rol);

	for(i=0; i < cur_rol; i++){
		cpd[i] = ipd[i+cur_id] - k_dt_date_to_count(pricing_date);
	}

	sum = pfixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv)* (*ZC).DF(spd/365.0,cpd[0]/365.);
	// printf("Fixed: %.10f %.10f\n",sum,(*UZC).DF(spd/365.0,cpd[0]/365.));
	for(i=1; i < cur_rol; i++){
		sum += pfixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)*(*ZC).DF(spd/365.0,cpd[i]/365.);
		// printf("Fixed: %.10f %.10f\n",pfixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)*(*UZC).DF(spd/365.0,cpd[i]/365.),(*UZC).DF(spd/360.0,cpd[i]/365.));
	}


	if (ipd) free(ipd);
	if (cpd) free(cpd);
	if (csd) free(csd);
	if (ced) free(ced);

	return(sum);

}

double k_swap_USDIRS_FLOATING_LEG(
								  date pricing_date,     // 가격산정일
								  date trade_date,       // trade date
								  date maturity,          // swap 총만기  
								  FRateType  frtype,     // Floating Rate Type (M3=3 month, M6=6 month)
								  double pfixed_rate,    // rate fixed previously
								  double margin,         // margin
								  int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
								  int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
								  k_cZeroCurve	*ZC // zero curve strucutre at pricing date
								  )
{
	vector<date> payment;  // vector for payment schedule
	int *ifd, *ipd;        // initial fixing schedule, payment schedule(거래 입력당시 확정된 스케쥴)
	int *cpd;              // current payment day-interval from pricing date (현재 남은 payment dates)
	int *csd, *ced;        // calculation start dates, end dates(거래 입력당시 확정되는 calculation period)
	int num_rol;           // initial number of rollovers(거래 입력 당시 확정되는 총 rollover 수)
	int cur_rol;           // current number of remaining rollovers(현 가격산정 시점에 남아있는 rollover 수)
	int cur_id;            // current index 
	int i;
	double sum;
	date spot_day;  // spot day 
	int spd;        // spot day interval from pricing day	
	int num_month;
	date effective_date;   // effective date

	if (frtype == M3)
		num_month = 3;
	else if (frtype == M6)
		num_month = 6;
	else if (frtype == M1)
		num_month = 1;


	num_rol = ((maturity.year-trade_date.year)*12+(maturity.month-trade_date.month))/num_month;

	effective_date = k_swap_getSpotDate(trade_date,NEW_YORK_LONDON);
	spot_day = k_swap_getSpotDate(pricing_date,NEW_YORK_LONDON);
	spd = k_dt_date_to_count(spot_day)-k_dt_date_to_count(pricing_date);

	ifd = (int *)malloc(sizeof(int)*num_rol);
	ipd = (int *)malloc(sizeof(int)*num_rol);
	csd = (int *)malloc(sizeof(int)*num_rol);
	ced = (int *)malloc(sizeof(int)*num_rol);

	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,NEW_YORK_LONDON);

	ipd[0] = k_dt_date_to_count(payment[0]);
	ifd[0] = k_dt_date_to_count(trade_date);
	csd[0] = k_dt_date_to_count(effective_date);
	ced[0] = ipd[0];

	i=0;

	for(i=1; i < num_rol; i++){
		ipd[i] = k_dt_date_to_count(payment[i]);
		ifd[i] = k_swap_day_shift(payment[i-1],LONDON,2,-1);  // 
		csd[i] = ipd[i-1];
		ced[i] = ipd[i];
	}


	for(i=0; i < num_rol; i++){

		if (k_dt_date_to_count(payment[i]) >= k_dt_date_to_count(pricing_date)) {
			cur_rol = num_rol-i;
			cur_id = i;
			break;
		}
	}


	cpd = (int *)malloc(sizeof(int)*cur_rol);

	for(i=0; i < cur_rol; i++){
		cpd[i] = ipd[i+cur_id] - k_dt_date_to_count(pricing_date);
	}


	double fra;
	vector<date> fraed;


	sum = (pfixed_rate+margin) * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv)* (*ZC).DF(spd/365.0,cpd[0]/365.);

	double dt1, dt2;
	date spotd;

	for(i=1; i < cur_rol; i++){

		spotd = k_swap_getSpotDate(k_dt_count_to_date(ifd[i+cur_id]),LONDON);
		fraed = k_swap_scheduleGenerator(spotd,1,num_month,BizDayConv,LONDON);

		dt1 = (k_dt_date_to_count(spotd)-k_dt_date_to_count(pricing_date))/365.0;
		dt2 = (k_dt_date_to_count(fraed[0])-k_dt_date_to_count(pricing_date))/365.0;

		fra = (*ZC).OR(dt1,dt2,k_dt_yearFrac(spotd,fraed[0],  DayCountConv));

		sum += (fra+margin)*k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)
			*(*ZC).DF(spd/365.0,cpd[i]/365.);
	}


	if (ifd) free(ifd);
	if (ipd) free(ipd);
	if (cpd) free(cpd);
	if (csd) free(csd);
	if (ced) free(ced);

	return(sum);

}

//	return USD PV for notional 1 USD
double k_swap_USDCRS_FLOATING_LEG(
								  date pricing_date,     // 가격산정일
								  date trade_date,       // trade date
								  date maturity,          // swap 총 만기
								  FRateType  frtype,     // Floating Rate Type (M3=3 month, M6=6 month)
								  double pfixed_rate,    // rate fixed previously
								  double margin,         // margin
								  int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
								  int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
								  calendar  pay_cal,
								  calendar  fix_cal,
								  k_cZeroCurve	*UZC // zero curve strucutre at pricing date
								  )
{
	vector<date> payment;  // vector for payment schedule
	int *ifd, *ipd;        // initial fixing schedule, payment schedule(거래 입력당시 확정된 스케쥴)
	int *cpd;              // current payment day-interval from pricing date (현재 남은 payment dates)
	int *csd, *ced;        // calculation start dates, end dates(거래 입력당시 확정되는 calculation period)
	int num_rol;           // initial number of rollovers(거래 입력 당시 확정되는 총 rollover 수)
	int cur_rol;           // current number of remaining rollovers(현 가격산정 시점에 남아있는 rollover 수)
	int cur_id;            // current index 
	int i;
	double sum=0.0;
	int spd;        // settlement day from pricing day	
	int num_month;
	date effective_date;   // effective date

	if (frtype == M3)
		num_month = 3;
	else if (frtype == M6)
		num_month = 6;
	else if (frtype == M1)
		num_month = 1;



	num_rol = ((maturity.year-trade_date.year)*12+(maturity.month-trade_date.month))/num_month;

	effective_date = k_dt_count_to_date(k_swap_day_shift(trade_date,pay_cal,2,1));
	spd = k_swap_day_shift(pricing_date,pay_cal,1,1)-k_dt_date_to_count(pricing_date);

	if ( k_dt_date_to_count(pricing_date) < k_dt_date_to_count(effective_date)){
		sum += (-1.0)*(*UZC).DF(spd/365.0, (k_dt_date_to_count(effective_date)-k_dt_date_to_count(pricing_date))/365.);

	}
	ifd = (int *)malloc(sizeof(int)*num_rol);
	ipd = (int *)malloc(sizeof(int)*num_rol);
	csd = (int *)malloc(sizeof(int)*num_rol);
	ced = (int *)malloc(sizeof(int)*num_rol);

	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,pay_cal);

	ipd[0] = k_dt_date_to_count(payment[0]);
	ifd[0] = k_dt_date_to_count(trade_date);
	csd[0] = k_dt_date_to_count(effective_date);
	ced[0] = ipd[0];
	for(i=1; i < num_rol; i++){
		ipd[i] = k_dt_date_to_count(payment[i]);
		ifd[i] = k_swap_day_shift(payment[i-1],fix_cal,2,-1);
		csd[i] = ipd[i-1];
		ced[i] = ipd[i];
	}


	for(i=0; i < num_rol; i++){
		if (k_dt_date_to_count(payment[i]) >= k_dt_date_to_count(pricing_date)) {
			cur_rol = num_rol-i;
			cur_id = i;
			break;
		}
	}


	cpd = (int *)malloc(sizeof(int)*cur_rol);
	for(i=0; i < cur_rol; i++){
		cpd[i] = ipd[i+cur_id] - k_dt_date_to_count(pricing_date);
	}


	double fra;
	vector<date> fraed;

	sum += (pfixed_rate+margin) * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv)* (*UZC).DF(spd/365.0,cpd[0]/365.);

	// printf("Floating: %.10f %.10f %.10f \n",pfixed_rate*100.,sum,(*UZC).DF(spd/365.0,cpd[0]/365.));

	double dt1, dt2;
	date spotd;

	for(i=1; i < cur_rol; i++){

		spotd = k_swap_getSpotDate(k_dt_count_to_date(ifd[i+cur_id]),fix_cal);
		fraed = k_swap_scheduleGenerator(spotd,1,num_month,BizDayConv,fix_cal);

		dt1 = (k_dt_date_to_count(spotd)-k_dt_date_to_count(pricing_date))/365.0;
		dt2 = (k_dt_date_to_count(fraed[0])-k_dt_date_to_count(pricing_date))/365.0;

		fra = (*UZC).OR(dt1,dt2,k_dt_yearFrac(spotd,fraed[0],  DayCountConv));
		sum += (fra+margin)*k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)
			*(*UZC).DF(spd/365.0,cpd[i]/365.);

		//printf("Floating: %.10f %.10f %.15f %d %d \n",fra*100.0,
		//	       fra*k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]), DayCountConv)*(*UZC).DF(spd/365.0,cpd[i]/365.),(*UZC).DF(spd/365.0,cpd[i]/365.),k_dt_date_to_count(fraed[0]),k_dt_date_to_count(spotd));
	}

	sum += (*UZC).DF(spd/365.0,cpd[cur_rol-1]/365.0);

	if (ifd) free(ifd);
	if (ipd) free(ipd);
	if (cpd) free(cpd);
	if (csd) free(csd);
	if (ced) free(ced);

	return(sum);

}


double k_swap_USDCRS_FIXED_LEG(
							   date pricing_date,     // 가격산정일
							   date trade_date,       // trade date
							   date maturity,         // 만기
							   PayType ptype,         // 이자지급주기 M=montly, Q=quarterly, S=semi-annually, A=annually 
							   double pfixed_rate,    // 확정고정금리
							   int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
							   int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
							   calendar pay_cal,
							   k_cZeroCurve	*UZC // zero curve strucutre at pricing date
							   )
{
	vector<date> payment;  // vector for payment schedule
	int *ipd;              // initial payment schedule(거래 입력당시 확정된 스케쥴)
	int *cpd;              // current payment day-interval from pricing date (현재 남은 payment dates)
	int *csd, *ced;        // calculation start dates, end dates(거래 입력당시 확정되는 calculation period)
	int num_rol;           // initial number of rollovers(거래 입력 당시 확정되는 총 rollover 수)
	int cur_rol;           // current number of remaining rollovers(현 가격산정 시점에 남아있는 rollover 수)
	int cur_id;            // current index 
	int i;
	double sum=0.0;
	int spd;        // settlement day from pricing day	
	int num_month;
	date effective_date;   // effective date

	if (ptype == M)
		num_month = 1;
	else if (ptype == Q)
		num_month = 3;
	else if (ptype == S)
		num_month = 6;
	else if (ptype == A)
		num_month = 12;


	num_rol = ((maturity.year-trade_date.year)*12+(maturity.month-trade_date.month))/num_month;
	effective_date = k_dt_count_to_date(k_swap_day_shift(trade_date,pay_cal,2,1));
	spd = k_swap_day_shift(pricing_date,pay_cal,1,1)-k_dt_date_to_count(pricing_date);

	if ( k_dt_date_to_count(pricing_date) < k_dt_date_to_count(effective_date)){
		sum += (-1.0)*(*UZC).DF(spd/365.0, (k_dt_date_to_count(effective_date)-k_dt_date_to_count(pricing_date))/365.);

	}

	ipd = (int *)malloc(sizeof(int)*num_rol);
	csd = (int *)malloc(sizeof(int)*num_rol);
	ced = (int *)malloc(sizeof(int)*num_rol);

	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,pay_cal);

	ipd[0] = k_dt_date_to_count(payment[0]);
	csd[0] = k_dt_date_to_count(effective_date);
	ced[0] = ipd[0];
	for(i=1; i < num_rol; i++){
		ipd[i] = k_dt_date_to_count(payment[i]);
		csd[i] = ipd[i-1];
		ced[i] = ipd[i];
		//	printf("%d %d %d \n",i,ipd[i],csd[i]);
	}

	for(i=0; i < num_rol; i++){

		if (k_dt_date_to_count(payment[i]) >= k_dt_date_to_count(pricing_date)) {
			cur_rol = num_rol-i;
			cur_id = i;
			break;
		}
	}

	cpd = (int *)malloc(sizeof(int)*cur_rol);

	for(i=0; i < cur_rol; i++){
		cpd[i] = ipd[i+cur_id] - k_dt_date_to_count(pricing_date);
	}

	sum += pfixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv)* (*UZC).DF(spd/365.0,cpd[0]/365.);
	// printf("0 %.10f %.15f %d %.15f \n",(*UZC).DF(spd/365.0,cpd[0]/365.),pfixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv)* (*UZC).DF(spd/365.0,cpd[0]/365.),ced[0]-csd[0],k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv));
	for(i=1; i < cur_rol; i++){
		sum += pfixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)*(*UZC).DF(spd/365.0,cpd[i]/365.);
		// printf("%d %.10f %.15f %d %.15f\n",i,(*UZC).DF(spd/365.0,cpd[i]/365.),pfixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)*(*UZC).DF(spd/365.0,cpd[i]/365.),ced[cur_id+i]-csd[cur_id+i],k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv));
	}

	sum += (*UZC).DF(spd/365.0,cpd[cur_rol-1]/365.0);

	if (ipd) free(ipd);

	if (cpd) free(cpd);
	if (csd) free(csd);
	if (ced) free(ced);

	return(sum);

}

//	return KRW PV for notioal 1 KRW
double k_swap_KRWCRS_FIXED_LEG(
							   date pricing_date,     // 가격산정일
							   date trade_date,       // trade date
							   date  maturity,         // swap 총 만기 
							   PayType ptype,         // 이자지급주기 M=montly, Q=quarterly, S=semi-annually, A=annually 
							   double fixed_rate,    // 확정고정금리
							   int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
							   int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
							   calendar  pay_cal,
							   k_cZeroCurve  *ZC  // zero curve strucutre at pricing date
							   )
{
	vector<date> payment;  // vector for payment schedule
	int *ipd;              // initial payment schedule(거래 입력당시 확정된 스케쥴)
	int *cpd;              // current payment day-interval from pricing date (현재 남은 payment dates)
	int *csd, *ced;        // calculation start dates, end dates(거래 입력당시 확정되는 calculation period)
	int num_rol;           // initial number of rollovers(거래 입력 당시 확정되는 총 rollover 수)
	int cur_rol;           // current number of remaining rollovers(현 가격산정 시점에 남아있는 rollover 수)
	int cur_id;            // current index 
	int i;
	double sum=0.0;
	int spd;        // settlement day from pricing day	
	int num_month;
	date effective_date;   // effective date


	if (ptype == M)
		num_month = 1;
	else if (ptype == Q)
		num_month = 3;
	else if (ptype == S)
		num_month = 6;
	else if (ptype == A)
		num_month = 12;


	num_rol = ((maturity.year-trade_date.year)*12+(maturity.month-trade_date.month))/num_month;
	//	printf("%d \n",num_rol);

	effective_date = k_dt_count_to_date(k_swap_day_shift(trade_date,pay_cal,2,1));
	spd = k_swap_day_shift(pricing_date,pay_cal,1,1)-k_dt_date_to_count(pricing_date);

	if ( k_dt_date_to_count(pricing_date) < k_dt_date_to_count(effective_date)){
		sum +=(-1.0)*(*ZC).DF(spd/365.0, (k_dt_date_to_count(effective_date)-k_dt_date_to_count(pricing_date))/365.);
	}

	ipd = (int *)malloc(sizeof(int)*num_rol);
	csd = (int *)malloc(sizeof(int)*num_rol);
	ced = (int *)malloc(sizeof(int)*num_rol);

	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,pay_cal);

	ipd[0] = k_dt_date_to_count(payment[0]);
	csd[0] = k_dt_date_to_count(effective_date);
	ced[0] = ipd[0];

	for(i=1; i < num_rol; i++){
		ipd[i] = k_dt_date_to_count(payment[i]);
		csd[i] = ipd[i-1];
		ced[i] = ipd[i];
		//	printf("%d %d %d \n",i,ipd[i],csd[i]);
	}

	for(i=0; i < num_rol; i++){
		if (k_dt_date_to_count(payment[i]) >= k_dt_date_to_count(pricing_date)) {
			cur_rol = num_rol-i;
			cur_id = i;
			break;
		}
	}

	cpd = (int *)malloc(sizeof(int)*cur_rol);
	for(i=0; i < cur_rol; i++){
		cpd[i] = ipd[i+cur_id] - k_dt_date_to_count(pricing_date);
	}

	sum += fixed_rate *  k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv)* (*ZC).DF(spd/365.0,cpd[0]/365.);
	for(i=1; i < cur_rol; i++){
		sum += fixed_rate * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)*(*ZC).DF(spd/365.0,cpd[i]/365.);
		//	   printf("%d %.10f \n",i,(*ZC).DF(spd/365.0,cpd[i]/365.));
	}

	sum += (*ZC).DF(spd/365.0,cpd[cur_rol-1]/365.0);


	if (ipd) free(ipd);
	if (cpd) free(cpd);
	if (csd) free(csd);
	if (ced) free(ced);

	return(sum);

}


double k_swap_KRWCRS_FLOATING_LEG(
								  date pricing_date,     // 가격산정일
								  date trade_date,       // trade date
								  date maturity,          
								  FRateType  frtype,     // Floating Rate Type (M3=3 month, M6=6 month)
								  double pfixed_rate,    // rate fixed previously
								  double margin,         // margin
								  int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
								  int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
								  calendar  pay_cal,     
								  calendar  fix_cal,
								  k_cZeroCurve  *FZC,   // zero curve structure for fixing 
								  k_cZeroCurve  *DZC   // zero curve structure for discounting 
								  )
{
	vector<date> payment;  // vector for payment schedule
	int *ifd, *ipd;        // initial fixing schedule, payment schedule(거래 입력당시 확정된 스케쥴)
	int *cpd;              // current payment day-interval from pricing date (현재 남은 payment dates)
	int *csd, *ced;        // calculation start dates, end dates(거래 입력당시 확정되는 calculation period)
	int num_rol;           // initial number of rollovers(거래 입력 당시 확정되는 총 rollover 수)
	int cur_rol;           // current number of remaining rollovers(현 가격산정 시점에 남아있는 rollover 수)
	int cur_id;            // current index 
	int i;
	double sum=0.0;
	int spd;        // settlement day from pricing day	
	int num_month;
	date effective_date;   // effective date

	if (frtype == M3)
		num_month = 3;
	else if (frtype == M6)
		num_month = 6;
	else if (frtype == M1)
		num_month = 1;



	num_rol = ((maturity.year-trade_date.year)*12+(maturity.month-trade_date.month))/num_month;

	//	printf("%d %d \n",  num_rol,num_month);

	effective_date = k_dt_count_to_date(k_swap_day_shift(trade_date,pay_cal,2,1));
	spd = k_swap_day_shift(pricing_date,pay_cal,1,1)-k_dt_date_to_count(pricing_date);

	if ( k_dt_date_to_count(pricing_date) < k_dt_date_to_count(effective_date)){
		sum +=(-1.0)*(*DZC).DF(spd/365.0, (k_dt_date_to_count(effective_date)-k_dt_date_to_count(pricing_date))/365.);
	}


	ifd = (int *)malloc(sizeof(int)*num_rol);
	ipd = (int *)malloc(sizeof(int)*num_rol);
	csd = (int *)malloc(sizeof(int)*num_rol);
	ced = (int *)malloc(sizeof(int)*num_rol);

	payment = k_swap_scheduleGenerator(effective_date,num_rol,num_month,BizDayConv,pay_cal);

	ipd[0] = k_dt_date_to_count(payment[0]);
	ifd[0] = k_dt_date_to_count(trade_date);
	csd[0] = k_dt_date_to_count(effective_date);
	ced[0] = ipd[0];
	i=0;
	//printf("%d %d %d %d %d \n",i,ipd[i],ifd[i],csd[i],ced[i]);
	for(i=1; i < num_rol; i++){
		ipd[i] = k_dt_date_to_count(payment[i]);
		ifd[i] = k_swap_day_shift(payment[i-1],fix_cal,2,-1);
		csd[i] = ipd[i-1];
		ced[i] = ipd[i];
		//printf("%d %d %d %d %d \n",i,ipd[i],ifd[i],csd[i],ced[i]);
	}


	for(i=0; i < num_rol; i++){

		if (k_dt_date_to_count(payment[i]) >= k_dt_date_to_count(pricing_date)) {
			cur_rol = num_rol-i;
			cur_id = i;
			break;
		}
	}


	cpd = (int *)malloc(sizeof(int)*cur_rol);

	for(i=0; i < cur_rol; i++){
		cpd[i] = ipd[i+cur_id] - k_dt_date_to_count(pricing_date);
	}


	double fra;
	vector<date> fraed;


	sum += (pfixed_rate+margin) * k_dt_yearFrac(k_dt_count_to_date(csd[cur_id]),k_dt_count_to_date(ced[cur_id]),  DayCountConv)* (*DZC).DF(spd/365.0,cpd[0]/365.);

	//printf("Floating: %.10f %.10f %.10f %d\n",pfixed_rate*100.,sum,(*FZC).DF(spd/365.0,cpd[0]/365.), k_dt_count_to_date(ced[cur_id])-k_dt_count_to_date(csd[cur_id]));


	double dt1, dt2;
	date spotd;

	for(i=1; i < cur_rol; i++){

		spotd = k_swap_getSpotDate(k_dt_count_to_date(ifd[i+cur_id]),fix_cal);
		fraed = k_swap_scheduleGenerator(spotd,1,num_month,BizDayConv,fix_cal);

		dt1 = (k_dt_date_to_count(spotd)-k_dt_date_to_count(pricing_date))/365.0;
		dt2 = (k_dt_date_to_count(fraed[0])-k_dt_date_to_count(pricing_date))/365.0;

		fra = (*FZC).OR(dt1,dt2,k_dt_yearFrac(spotd,fraed[0],  DayCountConv));
		sum += (fra+margin)*k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]),  DayCountConv)
			*(*DZC).DF(spd/365.0,cpd[i]/365.);

		//printf("Floating: %.10f %.10f %.15f %d %d \n",fra*100.0,
		//       fra*k_dt_yearFrac(k_dt_count_to_date(csd[cur_id+i]),k_dt_count_to_date(ced[cur_id+i]), DayCountConv)*(*FZC).DF(spd/365.0,cpd[i]/365.),(*FZC).DF(spd/365.0,cpd[i]/365.),k_dt_date_to_count(fraed[0]),k_dt_date_to_count(spotd));
	}

	sum += (*DZC).DF(spd/365.0,cpd[cur_rol-1]/365.0);
	//printf("%.10f \n",(*FZC).DF(spd/365.0,cpd[cur_rol-1]/365.0));

	if (ifd) free(ifd);
	if (ipd) free(ipd);
	if (cpd) free(cpd);
	if (csd) free(csd);
	if (ced) free(ced);

	return(sum);

}

