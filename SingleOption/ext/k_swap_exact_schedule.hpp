#ifndef kdb_swap_exact_schedule_hpp
#define kdb_swap_exact_schedule_hpp

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <iostream> //cout
using namespace std;

#include "kdb_BasicMath.hpp"
#include "kdb_LinearAlgebra.hpp"
#include "kdb_OptionFormula.hpp"
#include "kdb_probability.hpp"
#include "kdb_date.hpp"
#include "kdb_cZeroCurve.hpp"
#include "kdb_interpolation.hpp"
#include "kdb_miscellaneous.hpp"

/**********************************************************************/
/* Program Name      : kdb_swap_exact_schedule.hpp
/* version           : 0.0.0.1
/* author            : Jaemin Ahn, Kicheon Chang, Hangseop Cho
/* date              : 2009.12.26
/* Last modified by  : Jaemin Ahn
/* Last modified at  : 2010. 8. 24
/* Copyright         : KDB Quant team
/* Description       : functions involved in swap pricers with exact schedules
/* Related Doc. Name : none
/**********************************************************************/


enum FRateType{M1,M3,M6};    // Floating Rate Type M1 = 1 month, M3= 3 Month, M6 = 6 month
enum PayType{M=1,Q,S,A};  // M=monthly, Q=quarterly, S=semi-annually, A=annually


int k_swap_holiday_search(vector<long> hdays,int size,int key);

long k_swap_day_shift(date dt,calendar cal,int n, int direction);

/************************************************************************/
// shifting date to business date according to the business day convention
// BizDayConv: 1=following, 2=preceding, default=modified following
/************************************************************************/

date k_swap_dt_workday(date dt, int BizDayConv, calendar calendar_);

vector<date> k_swap_scheduleGenerator(date startDate,int numRollover,int schedule_term,int BizDayConv,calendar calendar_);

date k_swap_getSpotDate(date dt, calendar calendar_);

date k_swap_fix_date_shift(date dt, calendar calendar_);
// determine fixing day shifted back from payment date 

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
								  );

double k_swap_KRWIRS_FIXED_LEG(
							   date pricing_date,     // 가격산정일
							   date trade_date,       // trade date
							   date  maturity,         // swap 총 만기 unit:year
							   PayType ptype,         // 이자지급주기 M=montly, Q=quarterly, S=semi-annually, A=annually 
							   double pfixed_rate,    // 확정고정금리
							   int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
							   int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
							   k_cZeroCurve *ZC  // zero curve strucutre at pricing date
							   );

double k_swap_USDIRS_FIXED_LEG(
							   date pricing_date,     // 가격산정일
							   date trade_date,       // trade date
							   date  maturity,         // swap 총 만기 unit:year
							   PayType ptype,         // 이자지급주기 M=montly, Q=quarterly, S=semi-annually, A=annually 
							   double pfixed_rate,    // 확정고정금리
							   int BizDayConv,        // Business Day Convention: default(0)=modified following, 1=following, 2=preceding
							   int DayCountConv,      // Day Count Convention: default(0)=A/365, 1=30/360, 2=Actual/360, 3=Actual/Actual
							   k_cZeroCurve	*ZC // zero curve strucutre at pricing date
							   );

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
								  );

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
								  );

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
							   );

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
							   );

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
								  );

#endif