
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <new>
#include <algorithm>
#include <ios> 
#include "kdb_MonteCarlo.hpp"
#include "kdb_miscellaneous.hpp"
#include "kdb_probability.hpp"
#include "kdb_BasicMath.hpp"

using namespace std;

//======================================================================================
// Ran2 : pseudo-randorm number generator 
//            uniformly-distributed random numbers between 0 and 1 
//            code provided by Hangseob Cho
//            ref : Numerical Recipes (p.286)
//======================================================================================
double k_mc_ran2(int *idum)// when first called, it should be called by negative idum.
{
	const int IM1=2147483563, IM2=2147483399;
	const int IA1=40014,IA2=40692,IQ1=53668,IQ2=52774;
	const int IR1=12211,IR2=3791,NTAB=32,IMM1=IM1-1;
	const int NDIV=1+IMM1/NTAB;
	const double EPS=3.0e-16,RNMX=1.0-EPS,AM=1.0/double(IM1);
	static int idum2=123456789,iy =0;
	static int iv[NTAB];
	int j,k;
	double temp;

	if (*idum <= 0) {
		if (-(*idum)<1) *idum = 1;	// if (*idum)==0 or (*idum)== -(IM+1) then *idum = 1
		else *idum = -(*idum);
		idum2=*idum;
		for (j=NTAB+7;j>=0;j--) {	//load the shuffle table with random intergers(within [1,IM-1]) 
			k = (*idum)/IQ1;						
			*idum = IA1*(*idum-k*IQ1)-IR1*k;
			if(*idum<0) *idum += IM1;
			if (j<NTAB) iv[j] = *idum;
		}	
		iy = iv[0];
	}
	k = (*idum)/IQ1;
	*idum = IA1*(*idum-k*IQ1)-IR1*k;
	if (*idum<0) *idum += IM1;
	k = idum2/IQ2;
	idum2 = IA2*(idum2-k*IQ2)-IR2*k;
	if (idum2<0) idum2 += IM2;
	j = iy/NDIV;
	iy = iv[j]-idum2;
	iv[j] = *idum;
	if (iy <1) iy += IMM1;
	if ((temp=AM*iy)>RNMX) return RNMX;
	else return temp;
	
	//	!!caution!!: 
	//	when ran2 is first called it should be called by negative idum
	//	in order to make shuffle table.
	//	iy can get to be zero by successive call, so can't make new shuffle table
	//	on conditioon, iy == 0

	//	suspision:
	//	rand2의 idea는 기본적으로 period가 서로소인 ran0A와 ran0B를 더해서 period가 (period ran0A * period ran0A)가 되도록하는건데..
	//	문제는 이 코딩에서는 ran1A 와 ran0B를 더해서 구핸했다는것인데.. 
	//	ran1의 repeatition period를 정확히 모르는 상태인데 두개의 period가 서로소가 되는지 확신할 수 있는가?
}

//======================================================================================
//	NORMAL RANDOM NUMBER GENERATOR
//		added by Hangseob Cho at 2009-06-09
//======================================================================================
double k_mc_randn(int *idum)
{
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;

	if (*idum < 0) iset=0;
	if (iset==0) {
		do {
			v1 = 2.0*k_mc_ran2(idum)-1.0;	// ? idum or &idum
			v2 = 2.0*k_mc_ran2(idum)-1.0;	// ? idum or &idum
			rsq = v1*v1+v2*v2;
		} while (rsq>=1.0 || rsq==0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	}
	else {
		iset = 0;
		return gset;
	}
}

//======================================================================================
//	<<  SOBOL QUASI RANDOM NUMBER GENERATOR >>
// 
// INPUT: 
//   N         number of points  (cannot be greater than 2^32)
//   D         dimension      
//   dir_file  the input file containing direction numbers
//
// OUTPUT:
//
//     SOBOL[i][j] = the jth component of the ith point,
//   
//   with i indexed from 0 to N-1 and j indexed from 0 to D-1
//   
//   ***WARNING***  
//   PLEASE START FROM SOBOL[1][j] TO AVOID SOBOL[0][j]=0 
//   IN CASE YOU WANT TO GENERATE NORMALLY-DISTRIBUTED RANDOM NUMBERS
//	  FROM SOBOL NUMBERS USING INVERSE CUMULATIVE NORMAL FUNCTION
//
// a copy of the version by 
// Frances Y. Kuo
// Email: <f.kuo@unsw.edu.au>
//======================================================================================
double **k_mc_sobol(unsigned N, int D) //
{
  ifstream infile("kdb_sobol_direction_number.dat",ios::in);
  if (!infile) {
    cout << "Input file containing direction numbers cannot be found!\n";
    exit(1);
  }
  char buffer[1000];
  infile.getline(buffer,1000,'\n');
  infile.getline(buffer,1000,'\n');
  
  // L = max number of bits needed 
  unsigned L = (unsigned)ceil(log((double)N)/log(2.0)); 

  // C[i] = index from the right of the first zero bit of i
  unsigned *C = new unsigned [N];
  C[0] = 1;
  for (unsigned i=1;i<=N-1;i++) {
    C[i] = 1;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }
  
  // SOBOL[i][j] = the jth component of the ith point
  //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
  double **SOBOL = new double * [N];
  for (unsigned i=0;i<=N-1;i++) SOBOL[i] = new double [D];


  // ----- Compute the first dimension -----
  
  // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
  unsigned *V = new unsigned [L+1]; 
  for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1

 
  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  unsigned *X = new unsigned [N];
  X[0] = 0;
  SOBOL[0][0] = 0; 
  for (unsigned i=1;i<=N-1;i++) {
    X[i] = X[i-1] ^ V[C[i-1]];
    SOBOL[i][0] = (double)X[i]/pow(2.0,32); // *** the actual points
    //        ^ 0 for first dimension
  }

  // Clean up
  delete [] V;
  delete [] X;


  // ******* Compute the remaining dimensions ********
  for (int j=1;j<=D-1;j++) {    
    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    unsigned *m = new unsigned [s+1];
    for (unsigned i=1;i<=s;i++) infile >> m[i];
    
	// Compute direction numbers V[i], scaled by pow(2,32)
	// i = L for generating Sobol sequences but 
	// to check sobol properties, one may need to evaluate
	//further, depending on whether 2D is bigger than L. 
	//unsigned MAXDIR = L - 2*D >= 0 ? L : 2D;  
	unsigned *V = new unsigned [L+1]; 
	
    if (L <= s) {
      for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i); 
    }
    else {
      for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i); 
      for (unsigned i=s+1;i<=L;i++) {
		V[i] = V[i-s] ^ (V[i-s] >> s); 
		for (unsigned k=1;k<=s-1;k++) 
			V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
	  }
	} 

    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    unsigned *X = new unsigned [N];
    X[0] = 0 ;
	SOBOL[0][j] = 0; 
    for (unsigned i=1;i<=N-1;i++) {
      X[i] = X[i-1] ^ V[C[i-1]];
      SOBOL[i][j] = (double)X[i]/pow(2.0,32); // *** the actual points
      //        ^ j for dimension (j+1)
	}
	
    // Clean up
    delete [] m;
    delete [] V;
    delete [] X;
  } 
  //**************************************************

  delete [] C;

  
  return SOBOL;
}

//**************************//
// Calculate Sobol sequence //
//**************************//
double **k_mc_calculate_normalized_sobol_point(int Sobol_degree, int nb_stock, int nb_obs){

	int nb_SP	= power_HS( 2, Sobol_degree ) -1;		//	# Sobol Points
	int dim_SP 	= ( nb_obs ) * ( nb_stock );			//	dimension of Sobol Points

	double  **SP = k_mc_sobol( nb_SP + 1, dim_SP );
	double	**NSP = k_misc_new_double_array_HS( nb_SP, dim_SP );							

	if(SP && NSP){
		int i,j;
		for ( i = 0; i < nb_SP; i++ ) {
			for ( j = 0; j < dim_SP; j++ ) {
				NSP[i][j] = k_pr_ICN(SP[i+1][j]);
			}
		}
		k_misc_delete_array_HS( SP, nb_SP+1 );

	}
	if(SP && NSP)	return NSP;
	else			return NULL;
}
