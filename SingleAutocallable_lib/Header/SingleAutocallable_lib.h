#pragma once
void __stdcall SingleAutocallable_lib(double refprice, signed int expiryDate, int hitflag, int nb_autocall, signed int* auto_date,
	signed int vd, double spot,
	double* auto_strike, double* auto_coupon, double dummy_coupon, double kibarrier, double put_strike,
	double* Ivol, double* vol_term, double* vol_strike, int nb_vol_term, int nb_vol_strike,
	double* rfrate, double* rfrate_term, int nb_rfrate,
	double* divrate, double* divrate_term, int nb_divrate,
	double* result);

void __stdcall SingleAutocallable_discrete_dividend_lib(double refprice, signed int expiryDate, int hitflag, int nb_autocall, signed int* auto_date,
	signed int vd, double spot,
	double* auto_strike, double* auto_coupon, double dummy_coupon, double kibarrier, double put_strike,
	double* Ivol, double* vol_term, double* vol_strike, int nb_vol_term, int nb_vol_strike,
	double* rfrate, double* rfrate_term, int nb_rfrate,
	double* divrate, double* divrate_term, int nb_divrate,
	double* div_amount, signed int* div_exdate, int nb_div_amount,
	double* result);

void __stdcall SingleAutocallable_MC_lib(double refprice, signed int expiryDate, int hitflag, int nb_autocall, signed int* auto_date,
	signed int vd, double spot,
	double* auto_strike, double* auto_coupon, double dummy_coupon, double kibarrier, double put_strike,
	double* Ivol, double* vol_term, double* vol_strike, int nb_vol_term, int nb_vol_strike,
	double* rfrate, double* rfrate_term, int nb_rfrate,
	double* divrate, double* divrate_term, int nb_divrate,
	long numMc, double* result);

void __stdcall SingleAutocallable_MC_discrete_dividend_lib(double refprice, signed int expiryDate, int hitflag, int nb_autocall, signed int* auto_date,
	signed int vd, double spot,
	double* auto_strike, double* auto_coupon, double dummy_coupon, double kibarrier, double put_strike,
	double* Ivol, double* vol_term, double* vol_strike, int nb_vol_term, int nb_vol_strike,
	double* rfrate, double* rfrate_term, int nb_rfrate,
	double* divrate, double* divrate_term, int nb_divrate,
	double* div_amount, signed int* div_exdate, int nb_div_amount,
	long numMc, double* result);