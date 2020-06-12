#include "k_probability.hpp"

using namespace std;



double k_pr_erfccheb(double z)
// erf using Chebyshev coefficients
// z >= 0
{
	int j;
	double t,ty,tmp,d=0.0, dd = 0.0;

	t = 2./(2.+z);
	ty = 4.0*t-2.0;
	
	for(j=ncof-1; j>0; j--){
		tmp = d;
		d = ty*d - dd + cof[j];
		dd = tmp;
	}
	return t*exp(-z*z+0.5*(cof[0]+ty*d) - dd);
}

double k_pr_erf(double x)
// erf(x) error function
{
	if (x >= 0.0) 
		return 1.0-k_pr_erfccheb(x);
	else
		return k_pr_erfccheb(-x)-1.0;
}

double k_pr_erfc(double x)
// erfc(x) complementary error function
{
    if (x >= 0.0) 
		return k_pr_erfccheb(x);
	else
		return 2.0-k_pr_erfccheb(-x);
}

double k_pr_inverfc(double p)
// inverse of complementary error function
// 0 < p < 2
{
	double x, err, t, pp;
	int j;

	if (p >= 2.0) return -100.0;
	if (p <= 0.0) return 100.0;

	pp = (p < 1.0)? p : 2.0-p;
	t = sqrt(-2.*log(pp/2.0));
	x = -0.70711*((2.30753+t*0.27061)/(1.0+t*(0.99229+t*0.04481))-t);
	for(j=0; j < 2; j++){
		err = k_pr_erfc(x)-pp;
		x += err/(1.12837916709551257*exp(-x*x)-x*err);
	}
	return ( p < 1.0? x: -x);
}

double k_pr_inverf(double p)
// inverse of error function
// -1 < p < 1
{
	return( k_pr_inverfc(1.0-p));
}


double k_pr_LNpdf(double mu, double sig, double x)
// probability density function of Log-normal distribution
// x >= 0
{

	if (x == 0.0)
		return 0.0;
	return (0.398942280401432678/(sig*x))*exp(-0.5*(log(x)-mu)*(log(x)-mu)/sig/sig);
}

double k_pr_LNcdf(double mu, double sig, double x)
// cumulative distribution function of Log-normal distribution
// x >= 0
{
	if (x == 0.0)
		return 0.0;
	return 0.5*k_pr_erfc(-0.707106781186547524*(log(x)-mu)/sig);
}

double k_pr_LNinvcdf(double mu, double sig, double p)
// inverse cumulative distribution function of Log-normal distribution
// 0 < p < 1
{
	return exp(-1.41421356237309505*sig*k_pr_inverfc(2.0*p)+mu);
}


double k_pr_cdf_T_two_tailed(double t, double dof)
{
// returns the two-tailed cdf of student's t-distribution, prob(<|t|), t>0
	
	return 1.0-k_pr_ibeta(0.5*dof, 0.5, dof/(dof+t*t) );
}

double k_pr_cdf_F(double f, double dof1, double dof2)
{
// Returns cumulative distribution function of F distribution, prob(0<F<f)
	
	return k_pr_ibeta(0.5*dof1, 0.5*dof2, dof1*f/(dof2+dof1*f));
}


double k_pr_ibeta(double a, double b, double x)
{
// returns the incomplete beta function I_x(a,b) for positive a and b, and x between 0 and 1

	double bt;
	
	if (x == 0.0 || x == 1.0)		return x;
	if (a > SWITCH && b> SWITCH)	return k_pr_quad_ibeta(a,b,x);
	
	bt=exp( k_pr_ln_gamma(a+b) - k_pr_ln_gamma(a) - k_pr_ln_gamma(b) + a*log(x) + b*log(1.0-x) );

	if (x < (a+1.0)/(a+b+2.0))
		return bt*k_pr_cf_beta(a,b,x)/a;
	else
		return 1.0-bt*k_pr_cf_beta(b,a,1.0-x)/b;
}

double k_pr_cf_beta(const double a, const double b, const double x)
{
// Evaluates continued fraction for incomplete beta function by modified Lentz's method.

	const int MAXIT=10000;
	const double EPS=std::numeric_limits<double>::epsilon();
	const double FPMIN=std::numeric_limits<double>::min()/EPS;

	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++)
	{
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS) break;
	}
	// if (m > MAXIT) printf("a or b too big, or MAXIT too small in k_cf_beta");
	return h;
}

double k_pr_quad_ibeta(double a, double b, double x)
{
// Incomplete beta by quadrature. Returns I_x(a,b)
	const double y[18]={
		0.0021695375159141994,
		0.011413521097787704,
		0.027972308950302116,
		0.051727015600492421,
		0.082502225484340941,
		0.12007019910960293,
		0.16415283300752470,
		0.21442376986779355,
		0.27051082840644336,
		0.33199876341447887,
		0.39843234186401943,
		0.46931971407375483,
		0.54413605556657973,
		0.62232745288031077,
		0.70331500465597174,
		0.78649910768313447,
		0.87126389619061517,
		0.95698180152629142};
	const double w[18]={
		0.0055657196642445571,
		0.012915947284065419,
		0.020181515297735382,
		0.027298621498568734,
		0.034213810770299537,
		0.040875750923643261,
		0.047235083490265582,
		0.053244713977759692,
		0.058860144245324798,
		0.064039797355015485,
		0.068745323835736408,
		0.072941885005653087,
		0.076598410645870640,
		0.079687828912071670,
		0.082187266704339706,
		0.084078218979661945,
		0.085346685739338721,
		0.085983275670394821};

	int j;
	double xu, t, sum, ans;
	double a1=a-1.0, b1=b-1.0, mu=a/(a+b);
	double lnmu=log(mu), lnmuc=log(1.0-mu);
	t=sqrt(a*b/((a+b)*(a+b)*(a+b+1.0)));
	if( x > a/(a+b) )
	{
		if(x >=1.0) return 1.0;
		xu=min(1.0, max( mu + 10.0*t, x+ 5.0*t ));
	}
	else
	{
		if(x <= 0.0) return 0.0;
		xu=max(0.0, min(mu - 10.0*t, x-5.0*t));
	}
	sum=0.0;
	for(j=0;j<18;j++)
	{
		t=x + (xu-x)*y[j];
		sum += w[j]*exp(a1*(log(t)-lnmu) + b1*(log(1.0-t) - lnmuc));
	}
	ans=sum*(xu-x)*exp(a1*lnmu-k_pr_ln_gamma(a) + b1*lnmuc - k_pr_ln_gamma(b) + k_pr_ln_gamma(a+b));
	return ans>0.0? 1.0-ans : -ans;
}



double k_pr_ln_gamma(const double xx)
{
// Reurns the value ln[ Gamma(xx) ] for xx>0
	int j;
	double x,y,tmp,ser;
	static const double cof[14]={ 57.1562356658629235,
								-59.5979603554754912,
								 14.1360979747417471,
								 -0.491913816097620199,
								  0.339946499848118887e-4,
								  0.465236289270485756e-4,
								 -0.983744753048795646e-4,
								  0.158088703224912494e-3,
								 -0.210264441724104883e-3,
								  0.217439618115212643e-3,
								 -0.164318106536763890e-3,
								  0.844182239838527433e-4,
								 -0.261908384015814087e-4,
								  0.368991826595316234e-5};
	
	y=x=xx;
	tmp=x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser=0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}

void k_pr_gser(double *gamser, double a, double x, double *gln)
{
// Returns the incomplete gamma function P(a,x) evaluated by its series representation 
// as "gamser".  Also returns ln(Gamma(a)) as "gln"
//   x >= 0

	const int ITMAX=10000;
	const double EPS=std::numeric_limits<double>::epsilon();
	const double FPMIN=std::numeric_limits<double>::min()/EPS;

	int n;
	double sum,del,ap;

	*gln = k_pr_ln_gamma(a);
   
	if (x == 0.0){
		*gamser = 0.0;
		return;
	}
	else{
		ap = a;
		del = sum = 1.0/a;
		for(n=1; n=ITMAX; n++){
			ap++;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS){
				*gamser = sum * exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		printf("a too large, ITMAX too small in routine gser\n");
		return;
	}
	

}

void k_pr_gcf(double *gammcf, double a, double x, double *gln)
{
// Returns the incomplete gamma function Q(a,x) evaluated 
// by its continued fraction representation as "gamser".
// Also returns ln(Gamma(a)) as "gln"

	const int ITMAX=10000;
	const double EPS=std::numeric_limits<double>::epsilon();
	const double FPMIN=std::numeric_limits<double>::min()/EPS;


	int i;
	double an,b,c,d,del,h;

	*gln = k_pr_ln_gamma(a);
	b = x+1.0-a;
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d;
	
	for(i=1;i<=ITMAX;i++){
		an = -i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c = b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d = 1.0/d;
		del = d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) 
		printf("a too larte, ITMAX too small in gcf \n");
	*gammcf = exp(-x+a*log(x)-(*gln))*h;
}

double k_pr_igammap(double a, double x)
// Returns the incomplete gamma function P(a,x) 
// x >= 0, a > 0 
{

	 double gamser, gammcf, gln;
      
	 if (x < (a+1.0)){
		 k_pr_gser(&gamser,a,x,&gln);
		 return (gamser);
	 }
	 else{
		 k_pr_gcf(&gammcf,a,x,&gln);
		 return (1.0-gammcf);
	 }

}


double k_pr_igammaq(double a, double x)
// Returns the incomplete gamma function Q(a,x)=1-P(a,x) 
// x >= 0, a > 0 
{

	 double gamser, gammcf, gln;
      
	 if (x < (a+1.0)){
		 k_pr_gser(&gamser,a,x,&gln);
		 return (1.0-gamser);
	 }
	 else{
		 k_pr_gcf(&gammcf,a,x,&gln);
		 return (gammcf);
	 }

}


double k_pr_CP(double x, int k)
// Cumulative Poisson Probability Function P_x(< k) integer k >=1
// x : expected mean number 
{

	return(k_pr_igammaq((double)k,x));
}


double k_pr_chi_P(double chis, int dof)
// the probability that the observed chi-square for a correct model
//  should be less than a value chis
{
	return(k_pr_igammap(dof*0.5,chis*0.5));
}

double k_pr_chi_Q(double chis, int dof)
// the probability that the observed chi-square will exceed the value 
//   chis by chance even for a correct model.
{
    return(k_pr_igammaq(dof*0.5,chis*0.5));
}

double k_pr_Npdf(double mu, double sigma, double x)
// Probability density function for N(mu,sig)
{

	return( exp( - (x-mu)*(x-mu)*0.5/sigma/sigma)/sigma*ONE_OVER_SQRT_TWO_PI );
}



double k_pr_CN(double x)
//  Cummulative Univariate Normal Function
//  J. Hart(1968) Computer Approximations, Wiley. Algorithm 5666 for the error function
//  modified by G. West(2004)
//  : accurate to double precision throughout the real line
{
	  double xabs;
	  double cumnorm;
	  double build;
	  double exponential;
	  
	  xabs = fabs(x);
	  if (xabs > 37.0) {
	  	  cumnorm = 0.0;
	  }
	  else{
	  	  exponential = exp(-xabs*xabs*0.5);
	  	  if (xabs < 7.07106781186547){
	  	  	  build = (3.52624965998911e-02)*xabs + 0.700383064443688;
	  	  	  build = build*xabs + 6.37396220353165;
	  	  	  build = build*xabs + 33.912866078383;
	  	  	  build = build*xabs + 112.079291497871;
	  	  	  build = build*xabs + 221.213596169931;
	  	  	  build = build*xabs + 220.206867912376;
	  	  	  cumnorm = exponential * build;
	  	  	  build = (8.83883476483184e-02)*xabs + 1.75566716318264;
	  	  	  build = build*xabs + 16.064177579207;
	  	  	  build = build*xabs + 86.7807322029461;
	  	  	  build = build*xabs + 296.564248779674;
	  	  	  build = build*xabs + 637.333633378831;
	  	  	  build = build*xabs + 793.826512519948;
	  	  	  build = build*xabs + 440.413735824752;
	  	  	  cumnorm = cumnorm/build;
	    	}
	  	  else{
	  	  	  build = xabs + 0.65;
	  	  	  build = xabs + 4.0/build;
	  	  	  build = xabs + 3.0/build;
	  	  	  build = xabs + 2.0/build;
	  	  	  build = xabs + 1.0/build;
	  	  	  cumnorm = exponential/build/2.506628274631;
	  	  }
	  }
	  
	  if (x > 0)
	  	  return(1.0-cumnorm);
	  else
	  	  return(cumnorm);
	  	  
}

double k_pr_CBN(double x, double y, double rho)
// computing bivariate cumulative normal probabilities                        
// A. Genz (2004) 'Numerical computation of rectangular bivariate and trivariate normal and
//                 t probabilities', Statistics and Computing 14 251-260
// accurate to double precision
                   
{
	 int i, ISs, LG, NG;
	 double XX[11][4], W[11][4];
	 double h, k, hk, hs;
	 double BVN, ass, asr, sn;
	 double a, b, bs, c, d;
	 double xs, rs;
	 
	 W[1][1] = 0.17132449237917;
	 XX[1][1] = -0.932469514203152;
	 W[2][1] = 0.360761573048138;
	 XX[2][1] = -0.661209386466265;
	 W[3][1] = 0.46791393457269;
	 XX[3][1] = -0.238619186083197;
	
	 W[1][2] = 0.0471753363865118;
	 XX[1][2] = -0.981560634246719;
	 W[2][2] = 0.106939325995318;
	 XX[2][2] = -0.904117256370475;
	 W[3][2] = 0.160078328543346;
	 XX[3][2] = -0.769902674194305;
	 W[4][2] = 0.203167426723066;
	 XX[4][2] = -0.587317954286617;
	 W[5][2] = 0.233492536538355;
	 XX[5][2] = -0.36783149899818;
	 W[6][2] = 0.249147045813403;
	 XX[6][2] = -0.125233408511469;
	 
	 W[1][3] = 0.0176140071391521;
	 XX[1][3] = -0.993128599185095;
	 W[2][3] = 0.0406014298003869;
	 XX[2][3] = -0.963971927277914;
	 W[3][3] = 0.0626720483341091;
	 XX[3][3] = -0.912234428251326;
	 W[4][3] = 0.0832767415767048;
	 XX[4][3] = -0.839116971822219;
	 W[5][3] = 0.10193011981724;
	 XX[5][3] = -0.746331906460151;
	 W[6][3] = 0.118194531961518;
	 XX[6][3] = -0.636053680726515;
	 W[7][3] = 0.131688638449177;
	 XX[7][3] = -0.510867001950827;
	 W[8][3] = 0.142096109318382;
	 XX[8][3] = -0.37370608871542;
	 W[9][3] = 0.149172986472604;
	 XX[9][3] = -0.227785851141645;
	 W[10][3] = 0.152753387130726;
	 XX[10][3] = -0.0765265211334973;
	 
	 if (fabs(rho) < 0.3){
	 	   NG = 1;
	 	   LG = 3;
	 }
	 else if (fabs(rho) < 0.75){
	 	       NG = 2;
	 	       LG = 6;
	 	    }
	 	    else{
	 	    	 NG = 3;
	 	    	 LG = 10;
	 	    }
	 
	 h = -x;
	 k = -y;
	 hk = h*k;
	 BVN = 0.0;
	 	    
	 if (fabs(rho) < 0.925){
	 	   if (fabs(rho) > 0.0 ){
	 	   	   hs = (h*h+k*k)*0.5;
	 	   	   asr = asin(rho);
	 	   	   for(i=1; i <= LG; i++){
	 	   	   	   for(ISs = -1; ISs <=-1; ISs += 2){
	 	   	   	   	   sn = sin(asr*(ISs*XX[i][NG]+1.0)*0.5);
	 	   	   	   	   BVN += W[i][NG]*exp((sn*hk-hs)/(1.0-sn*sn));
	 	   	   	   }
	 	   	   }
	 	   	   BVN *= asr/(4.0*PI);
	 	   }
	 }
	 else{
	 	   if (rho < 0){
	 	   	   k = -k;
	 	   	   hk = -hk;
	 	   }
	 	   if (fabs(rho) < 1.0){
	 	   	   ass = (1.0-rho)*(1.0+rho);
	 	   	   a = sqrt(ass);
	 	   	   bs = (h-k)*(h-k);
	 	   	   c = (4.0-hk)/8.0;
	 	   	   d = (12.0-hk)/16.0;
	 	   	   asr = (-bs/ass + hk)*0.5;
	 	       if (asr > -100.0){
	 	       	   BVN = a*exp(asr)*(1.0-c*(bs-ass)*(1.0-d*bs/5.0)/3.0 + c*d*ass*ass/5.0);
	 	       }
	 	       if (-hk < 100.0){
	 	       	   b = sqrt(bs);
	 	       	   BVN -= exp(-hk*0.5)*sqrt(2.0*PI)*k_pr_CN(-b/a)*b*(1.0-c*bs*(1.0-d*bs/5.0)/3.0);
	 	       }
	 	       a *= 0.5;
	 	       for(i=1; i <= LG; i++){
	 	       	   for(ISs = -1; ISs <= 1; ISs += 2){
	 	       	   	   xs = (a*(ISs*XX[i][NG]+1.0))*(a*(ISs*XX[i][NG]+1.0));
	 	       	   	   rs = sqrt(1.0-xs);
	 	       	   	   asr = -(bs/xs + hk)*0.5;
	 	       	   	   if (asr > -100.0){
	 	       	   	   	   BVN += a*W[i][NG]*exp(asr)*(exp(-hk*(1.0-rs)/(2.0*(1.0+rs)))/rs
	 	       	   	   	          - (1.0 + c*xs*(1.0+d*xs)));
	 	       	   	   }
	 	       	   }
	 	       }
	 	       BVN = -BVN/(2.0*PI);
	 	   }
	 	   
	 	   if (rho > 0.0)
	 	   	   BVN += k_pr_CN(-max(h,k));
	 	   else{
	 	   	   BVN = -BVN;
	 	   	   if (k > h) 
	 	   	      BVN += k_pr_CN(k)-k_pr_CN(h);
	 	   }	      
	 }
	 return(BVN);
}



//
double k_pr_ICN(double prob)  //  0 < prob < 1
// The Inverse cumulative normal distribution function
// by Peter J. Acklam, University of Oslo, Statistics Division.
//
// URL: http://www.math.uio.no/~jacklam/notes/invnorm
{   

  // Coefficients for the rational approximation.
  const double
    a1 = -3.969683028665376e+01,
    a2 =  2.209460984245205e+02,
    a3 = -2.759285104469687e+02,
    a4 =  1.383577518672690e+02,
    a5 = -3.066479806614716e+01,
    a6 =  2.506628277459239e+00;
    
  const double
    b1 = -5.447609879822406e+01,
    b2 =  1.615858368580409e+02,
    b3 = -1.556989798598866e+02,
    b4 =  6.680131188771972e+01,
    b5 = -1.328068155288572e+01;
    
  const double
    c1 = -7.784894002430293e-03,
    c2 = -3.223964580411365e-01,
    c3 = -2.400758277161838e+00,
    c4 = -2.549732539343734e+00,
    c5 =  4.374664141464968e+00,
    c6 =  2.938163982698783e+00;
    
  const double
    d1 =  7.784695709041462e-03,
    d2 =  3.224671290700398e-01,
    d3 =  2.445134137142996e+00,
    d4 =  3.754408661907416e+00;
    

  // Limits of the approximation region.
  const double
    u_low   = 0.02425,
    u_high  = 1.0 - u_low;

  

  register double z, r;

  // Rational approximation for the lower region. ( 0 < prob < u_low )
  if( prob < u_low ){
    z = sqrt(-2.0*log(prob));
    z = (((((c1*z+c2)*z+c3)*z+c4)*z+c5)*z+c6) / ((((d1*z+d2)*z+d3)*z+d4)*z+1.0);
  }
  // Rational approximation for the central region. ( u_low <= prob <= u_high )
  else if( prob <= u_high ){
    z = prob - 0.5;
    r = z*z;
    z = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*z / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0);
  }
  // Rational approximation for the upper region. ( u_high < prob < 1 )
  else {
    z = sqrt(-2.0*log(1.0-prob));
    z = -(((((c1*z+c2)*z+c3)*z+c4)*z+c5)*z+c6) /  ((((d1*z+d2)*z+d3)*z+d4)*z+1.0);
  }

#define REFINE_INVERSECUMULATIVENORMAL_TO_FULL_MACHINE_PRECISION_USING_HALLEYS_METHOD

#ifdef REFINE_INVERSECUMULATIVENORMAL_TO_FULL_MACHINE_PRECISION_USING_HALLEYS_METHOD
  // The relative error of the approximation has absolute value less
  // than 1.15e-9.  One iteration of Halley's rational method (third
  // order) gives full machine precision.

  r = (k_pr_CN(z) - prob) * SQRT_TWO_PI * exp( 0.5 * z * z );	//	f(z)/df(z)
  z -= r/(1+0.5*z*r);							//	Halley's method
#endif

  return z;
}
