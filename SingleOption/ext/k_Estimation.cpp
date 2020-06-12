#include "kdb_Estimation.hpp"
#include "kdb_probability.hpp"
/**********************************************************************/
/* Program Name      : regress
/* version           : 0.1.0.0
/* author            : Jihyun Lee
/* date              : 2009. 1. 7
/* Modified by       : Jihyun Lee
/* Modified at       : 
/* Copyright         : KDB Quant team
/* Description       : Ordinary Least Square Regression.
/* Related Doc. Name : EViews 6 Users Guide II.pdf pp.11-16
/**********************************************************************/
/* PURPOSE:
*		return the OLS estimates, b_ols for the regression equation y=x*b+epsilon
* 
*  INPUTS:
*		y:  A N by 1 matrix, dependent variable
*		x:	A N by K matrix, explanation variables, x[i][0]=1.0 for all i.	
* 
*  OUTPUTS:
*		result:	A K by 4 matrix, Results for the OLS regression
*				column1: coefficient estimates, bOLS
*				column2: standard errors of the estimates
*				column3: t-statistics used to test the hypothesis that a coefficiet is equal to zero
*				column4: Prob(>|t|)
*		R2:		the R-squared statistic
*		adjR2:	the adjusted R2
*		F:		F-statistic used to test the hypothesis that all of the slope coefficients in a regression are zero.
*		probF:	Prob(>F)
*
*  USAGE:
*		void main()
*		{
*			int i;
*			double R2,adjR2, F, probF;	// diagnostic statistics
*			matrix<double> y(6,1);		// dependent variable
*			matrix<double> x(6,2);		// explanation variables
*
*			for(i=0;i<6;i++)	x[i][0]=1.0;
*			x[0][1]=1.0;	x[1][1]=2.3;	x[2][1]=3.1;	x[3][1]=4.8;	x[4][1]=5.6;	x[5][1]=6.3;
*			y[0][0]=2.6;	y[1][0]=2.8;	y[2][0]=3.1;	y[3][0]=4.7;	y[4][0]=5.1;	y[5][0]=5.3;
*			
*			matrix<double> result=k_est_regress(y,x,&R2, &adjR2, &F, &probF);
*			cout << "Estimates\tStd.Err\tt-Statistic\tProb(>|t|)" << endl;
*			k_mx_cout(result);
*
*			cout << "R2=\t" << R2 <<endl;
*			cout << "adjR2=\t" << adjR2 << endl;
*			cout << "F=\t" << F << endl;
*			cout << "probF=\t" << probF << endl;
*		}
*		*****************************************************************
*		Estimates   Std.Err		t-Statistic Prob(>|t|)
*		1.68422		0.290559    5.79649		0.00440436
*		0.584184    0.0678631   8.60828		0.0010009
*		R2=     0.948785
*		adjR2=  0.935982
*		F=      74.1024
*		probF=  0.0010009
*
/**********************************************************************/


matrix<double> k_est_regress(matrix<double> &y, matrix<double> &x, double *R2, double *adjR2, double *F, double *probF )
{
	int i;
	int N=x.numrows();	// number of observation
	int K=x.numcols();	// number of explanation variables+1

	matrix<double> bols	= k_mx_inverse(k_mx_transpose(x)*x)*k_mx_transpose(x)*y;	// ols estimates
	matrix<double> e	= y - x*bols;										// estimation errors
	matrix<double> ssr	= k_mx_transpose(e)*e;								// sum of squared residuals
	
	double s2 = ssr[0][0]/(N-K);											// estimator of unknown variance of errors
	matrix<double> varbols = s2*k_mx_inverse(k_mx_transpose(x)*x);				// variance-covariance matrix of bols

	matrix<double> Tstats(K,1);												// t-statistics
	matrix<double> probT(K,1);												// prob(>|t|)
	for(i=0;i<K;i++)
	{
		Tstats[i][0]= bols[i][0]/sqrt(varbols[i][i]);
		probT[i][0]	= 1.0 - k_pr_cdf_T_two_tailed(abs(Tstats[i][0]),N-K);
	}
	
	double ybar=0.0;														// ybar: mean of y
	for(i=0;i<N;i++)	ybar += y[i][0];
	ybar/=(double)N;
	matrix<double> mxybar	= ybar*k_mx_unit_matrix(N,1);
	matrix<double> meandevy	= y - mxybar;									// mean deviated form of y
	matrix<double> yprimey	= k_mx_transpose(meandevy)*meandevy;				

	*R2		= 1.0 - ssr[0][0]/yprimey[0][0];
	*adjR2	= 1.0 - (1.0-*R2)*(N-1)/(N-K);
	*F		= *R2 * (N-K)/(1-*R2)/(K-1);
	*probF	= 1 - k_pr_cdf_F(*F, K-1, N-K);

	matrix<double> result(K,4);
	for(i=0;i<K;i++)
	{
		result[i][0]=bols[i][0];
		result[i][1]=sqrt(varbols[i][i]);
		result[i][2]=Tstats[i][0];
		result[i][3]=probT[i][0];
	}

	return result;
}



