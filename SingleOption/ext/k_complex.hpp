#ifndef KDB_COMPLEX_H
#define KDB_COMPLEX_H

#include <math.h>
#include <algorithm>

using namespace std;

/**********************************************************************/
/* Program Name      : COMPLEX C library
/* version           : 0.1.0.0
/* author            : Jaemin Ahn
/* date              : 2009. 9. 15
/* Modified by       : 
/* Modified at       : 
/* Copyright         : KDB Quant team
/* Description       : complex type, operations etc
/* Related Doc. Name : none
/**********************************************************************/



typedef struct{   /* define the data type for complex numbers */
  double re;
  double im;
} k_complex;


double k_com_abs(k_complex z);   // absolute value of a complex number z

k_complex k_com_rmul(double r, k_complex z);  // (real number r) x (complex number z)

k_complex k_com_exp(k_complex z);            // complex exponential exp(z)

k_complex k_com_mul(k_complex z1,k_complex z2);  // multiplication of two complex numbers z1 & z2 

k_complex k_com_add(k_complex z1,k_complex z2);  // addition z1 + z2 of two complex numbers z1 & z2 

k_complex k_com_sub(k_complex z1,k_complex z2);  // subtraction z1 - z2 of two complex numbers z1 & z2 

k_complex k_com_div(k_complex z1,k_complex z2);  // division z1/z2  of two complex numbers z1 & z2 

k_complex k_com_sqrt(k_complex z);  // square root sqrt(z) of a complex number z 

k_complex k_com_inv(k_complex z);   // inverse 1/z of a complex number z

k_complex k_com_log(k_complex z);   // logarithm log(z) of a complex number 

#endif