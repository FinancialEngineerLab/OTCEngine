#ifndef KDB_DATE_H
#define KDB_DATE_H

#include <cmath>
#include <ctime>
#include <vector>
using namespace std;

/**********************************************************************/
/* Program Name      : date library
/* version           : 1.0.0.1
/* author            : kicheon chang
/* date              : 2008.12.24
/* Last Modified by  : Great Quant Jaemin Ahn
/* Last Modified at  : 2011.01.17
/* Copyright         : KDB Quant team
/* Description       : functions related to dates, days
/* Related Doc. Name : none
/**********************************************************************/

#define TRUE	1
#define FALSE	0

enum {JAN=1,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC};
enum calendar{SEOUL=1, NEW_YORK, LONDON, TARGET, TOKYO, NEW_YORK_LONDON, 
SEOUL_NEW_YORK, TARGET_NEW_YORK, TOKYO_NEW_YORK };

/************************************************************************/
// date class(structure)
//
/************************************************************************/
class date
{
public:
	int year, month, day;
	void operator=(int iYYYYMMDD)
	{
		this->year=int(iYYYYMMDD/10000);
		this->month=int((iYYYYMMDD-year*10000)/100);
		this->day=(iYYYYMMDD-year*10000-month*100);	
	};
};
//Operator
date operator+(const date &a, long b);
date operator-(const date &a, long b);
int operator-(const date &a, const date &b);

//Converting function
long			k_dt_date_to_count(date dt);

date			k_dt_count_to_date(long iCount);
date			k_dt_yyyymmdd_to_date(int YYYYMMDD);
int				k_dt_date_to_yyyymmdd(date dt);

//Information function
long			k_dt_week(long base_date_count);
int				k_dt_weekDate(int iYYYYMMDD);
int				k_dt_weekDate(date dt);
date			k_dt_Today();
double			k_dt_yearFrac(date dt1, date dt2, int DayCountConv);
double			k_dt_yearFrac(int dt1,int dt2, int DayCountConv);

int				k_dt_numDays(date dt1, date dt2, int DayCountConv);
date			k_dt_getSpotDate(date dt, calendar calendar_);
void			k_dt_getHolidays(vector<long> &hdays,calendar calendar_ );
int				k_dt_getLastDayOfMonth(date dt);
vector<date>	k_dt_scheduleGenerator(date startDate,int numRollover,int schedule_term,int BizDayConv,calendar calendar_);
int				k_dt_is_leap_year(int year);
int				k_dt_is_bizDay(date dt, calendar calendar_);

//Date-Shifting function
date			k_dt_month_adder(date dt, long add_month);
date			k_dt_workday(date dt,int BizDayConv, calendar calendar_);
long			k_dt_day_shift(date dt,calendar cal,int n, int direction);
vector<int>     k_dt_vday_shift(vector<int> dts,calendar cal, int n, int direction);

//internal use only
long			k_dt_workday(long base_date, long direction); //check Sunday and Saturday
int				k_dt_holiday_search(vector<long> hdays,int size,int key);
date			k_dt_workday_intrim(date dt,int direction,calendar calendar_);
void			k_dt_getHolidays_SEOUL  (vector<long> &hdays);
void			k_dt_getHolidays_NEW_YORK  (vector<long> &hdays);
void			k_dt_getHolidays_LONDON  (vector<long> &hdays);
void			k_dt_getHolidays_TARGET  (vector<long> &hdays);
void			k_dt_getHolidays_TOKYO  (vector<long> &hdays);
void			k_dt_getHolidays_NEW_YORK_LONDON  (vector<long> &hdays);
void			k_dt_getHolidays_SEOUL_NEW_YORK  (vector<long> &hdays);
void			k_dt_getHolidays_TARGET_NEW_YORK  (vector<long> &hdays);
void			k_dt_getHolidays_TOKYO_NEW_YORK  (vector<long> &hdays);
#endif