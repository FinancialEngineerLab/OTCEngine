#ifndef KDB_ESTIMATION_H
#define KDB_ESTIMATION_H

#include <vector>
#include <math.h>
#include "kdb_LinearAlgebra.hpp"

using namespace std;

/**********************************************************************/
/* Program Name      : Estimation library
/* version           : 0.1.0.0
/* author            : Jihyun Lee
/* date              : 2009. 1. 7
/* Modified by       : Jihyun Lee
/* Modified at       : 
/* Copyright         : KDB Quant team
/* Description       : Ordinary Least Square Regression
/* Related Doc. Name : none
/**********************************************************************/

matrix<double> k_est_regress(matrix<double> &y, matrix<double> &x, double *R2, double *adjR2, double *F, double *probF );

#endif