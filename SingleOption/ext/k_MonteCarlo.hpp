#ifndef KDB_MONTE_CARLO_HPP
#define KDB_MONTE_CARLO_HPP

/*************************************************************************
* Program Name : kdb_MonteCarlo.hpp
* Version : 0.0.0.0
* AuthoRate : Hwang-hyun Kwon
* Date : Jan. 15, 2009
* Modified by : Hwang-hyun Kwon
* Modified at : Jan. 15, 2009
* Copyright : KDB Quantitative Analysis Team 
* Description : A collection of pseudo and quasi random number generators
**************************************************************************/

double k_mc_ran2(int *iseed); 

double k_mc_randn(int *idum);

double **k_mc_sobol(unsigned uNumSimul , int iDimension);

double **k_mc_calculate_normalized_sobol_point(int Sobol_degree, int nb_stock, int nb_RSV);
						 
#endif
