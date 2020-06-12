#include "k_interpolation.hpp"
#include <cmath>
#include <stdexcept>

/************************************************************************/
// spline interpolation
// iDim; 1=linear, 3=cubic
/************************************************************************/
double k_intp_spline(vector<double>x, vector<double> y, double x_star, int iDim)
{
	switch(iDim)
	{	
	case 1: 
		return k_intp_linear_spline(x,y,x_star);
	case 3:
		return k_intp_cubic_spline(x,y,x_star);
	default:
		return k_intp_cubic_spline(x,y,x_star);
	}
}

/************************************************************************/
// 2 dimensional linear interpolation
//
/************************************************************************/
double k_intp_bi_linear(double x1, double x2, double y1, double y2, 
					  double z11, double z12, double z21, double z22, double x_star,double y_star)
{
	double z1121, z1222, z_star;
	z1121=k_intp_linear(x1,x2,z11,z21,x_star);
	z1222=k_intp_linear(x1,x2,z12,z22,x_star);
	z_star=k_intp_linear(y1,y2,z1121,z1222,y_star);
	return z_star;
}

/************************************************************************/
// linear spline interpolation
//
/************************************************************************/
double k_intp_linear_spline(vector<double> x, vector<double> y, double x_star)
{
	int i;
	double result;
	if(x.size()!=y.size())
		return 0;
	for(i=0;i<(int)x.size();i++)
	{
		if(x_star<=x[i])
			break;
	}
	result=k_intp_linear(x[i],x[i+1],y[i],y[i+1],x_star);

	return result;
}


/************************************************************************/
// polynomial interpolation
//
/************************************************************************/
double k_intp_polynomial(vector<double> x, vector<double> y, double x_star)
{
	double dSum, dProduct;
	int i,j;
	if(x.size()!=y.size())
		return 0;
	
	dSum=0;
	for(i=0;i<(int)x.size();i++)
	{
		dProduct=y[i];
		for(j=0;j<(int)x.size();j++)
		{
			if(i!=j)
				dProduct=dProduct*(x_star-x[j])/(x[i]-x[j]);
		}
		dSum=dSum+dProduct;
	}
	return dSum;
}

/************************************************************************/
// cubic spline interpolation
//
/************************************************************************/
double k_intp_cubic_spline(vector<double> x, vector<double> y, double x_star)
{
	int n;
	n=(int)x.size();
	if(x.size()!=y.size())
		return 0;

	vector<double> h(n);
	vector<double> a(n);
	vector<double> b(n);
	vector<double> c(n);
	vector<double> r(n);
	vector<double> y2(n);
	vector<double> coef_a(n);
	vector<double> coef_b(n);
	vector<double> coef_c(n);
	vector<double> coef_d(n);


	int i;
	double temp;
	double value;
	int m=n-1;

	for(i=0;i<m;i++)
	{
		h[i]=x[i+1]-x[i];
	}
	for(i=1;i<m;i++)
	{
		a[i]=h[i-1];
		b[i]=2.0*(h[i-1]+h[i]);
		c[i]=h[i];
	}
	for(i=1;i<m;i++)
	{
		r[i]=6.0*((y[i+1]-y[i])/h[i]+(y[i-1]-y[i])/h[i-1]);
	}

	y2[0]=0;
	y2[m]=0;
	a[1]=0;
	c[m-1]=0;

	for(i=1;i<m;i++)
	{
		temp=a[i+1]/b[i];
		b[i+1]=b[i+1]-temp*c[i];
		r[i+1]=r[i+1]-temp*r[i];
	}

	y2[m-1]=r[m-1]/b[m-1];
	for(i=m-2;i>=1;i--)
	{
		y2[i]=(r[i]-c[i]*y2[i+1])/b[i];
	}

	for(i=1;i<=m;i++)
	{
		coef_a[i]=y2[i-1]/(6.0*h[i-1]);
		coef_b[i]=y2[i]/(6.0*h[i-1]);
		coef_c[i]=y[i-1]/h[i-1]-y2[i-1]*h[i-1]/6.0;
		coef_d[i]=y[i]/h[i-1]-y2[i]*h[i-1]/6.0;
	}

	for(i=1;i<=m;i++)
	{
		if (x_star>=x[i-1] && x_star<=x[i])
		{value=coef_a[i]*pow((x[i]-x_star),3)
		+coef_b[i]*pow((x_star-x[i-1]),3)
		+coef_c[i]*(x[i]-x_star)
		+coef_d[i]*(x_star-x[i-1]);
		break;		
		}
		value=0.0;
	}
	return value;
}

/************************************************************************/
// linear interpolation between 2 points
//
/************************************************************************/
double k_intp_linear(double x1, double x2, double y1, double y2, double x_star)
{
	double value;
	double dx,dy;
	dx=(x2-x1);
	dy=(y2-y1);
	value=y1+dy/dx*(x_star-x1);
	return value;
}

/************************************************************************/
// Extraploate
// ext_method; 0=linear, 1=flat, 2=function_extension
/************************************************************************/
double k_intp_extrapolate(vector<double> x, vector<double> y, double x_star, int ext_method)
{

	int n=(int)x.size();
	switch(ext_method)
	{
	case 0: //linear extrapolation
		if(x_star<=x[0])
			return k_intp_linear(x[0],x[1],y[0],y[1],x_star);
		else if(x_star>=x[n-1])
			return k_intp_linear(x[n-2],x[n-1],y[n-2],y[n-1],x_star);
		else
			return k_intp_linear_spline(x,y,x_star); //lnterpolation
	case 1:
		if(x_star<=x[0])
			return y[0];
		else if(x_star>=x[n-1])
			return y[n-1];
		else
			return k_intp_linear_spline(x,y,x_star); //lnterpolation
	case 2:
/*		if(x_star<=x[0])
			return y[0];
		else if(x_star>=x[n-1])
			return y[n-1];
		else
			return k_intp_linear_spline(x,y,x_star); //lnterpolation
*/
		return 0.0;
	default:
		if(x_star<=x[0])
			return k_intp_linear(x[0],x[1],y[0],y[1],x_star);
		else if(x_star>=x[n-1])
			return k_intp_linear(x[n-2],x[n-1],y[n-2],y[n-1],x_star);
		else
			return k_intp_linear_spline(x,y,x_star); //lnterpolation
		}
}

/////////////////////////
//	added by hangseob 090331

double	k_intp_bilin_interp( 
						double xm, double ym,

						double xd, double xu,
						double yd, double yu,

						//	V(x,y)'s
						double Vdd, double Vdu,		
						double Vud, double Vuu		
						)
{	
	return (
		Vdd*(xu-xm)*(yu-ym) +
		Vdu*(xu-xm)*(ym-yd) +
		Vud*(xm-xd)*(yu-ym) +
		Vuu*(xm-xd)*(ym-yd))
		/(xu-xd)/(yu-yd);
}

double	k_intp_bilin_interp(
						double	xMin,	double	xMax, 	int 	nx,
						double	yMin,	double	yMax, 	int 	ny,
						double	**fxy,	double	x,	  	double 	y		)
{
	int		LI_x;	//	lower index of x
	int		LI_y;	//	lower index of y
	double	dx, dy;

	dx 		= (xMax - xMin)/nx;
	dy 		= (yMax - yMin)/ny;
	LI_x	= (int)floor( (x - xMin)/dx );
	LI_y	= (int)floor( (y - yMin)/dy );	

	if ( LI_x <= 0    ) LI_x = 0;
	if ( LI_x >= nx-1 ) LI_x = nx-1;
	if ( LI_y <= 0    ) LI_y = 0;
	if ( LI_y >= ny-1 ) LI_y = ny-1;

	return	k_intp_bilin_interp(	x,					y,
		xMin + dx*LI_x, 	xMin + dx*(LI_x+1),	
		yMin + dy*LI_y,		yMin + dy*(LI_y+1),
		fxy[LI_x  ][LI_y],	fxy[LI_x  ][LI_y+1],
		fxy[LI_x+1][LI_y],	fxy[LI_x+1][LI_y+1] 	);
}	

double k_intp_linear( double x_min, double x_max, int nx, double *fx, double x_star )
//========================================================================================
//	added by Hangseob at 090603
//	xmin = x[0], x[1], ... , x[nx] = xmax have function values fx[0], fx[1], ... , fx[nx]
//	x[i]'s are increasing and evenly spaced
//	I.E. 
//	fx[0]   = f(x_min)
//	fx[nx]	= f(x_max)
//	fx[k]   = f(x_min + (x_max - x_min)/nx*k) 
//
//	The function returns f(x_star) by linear interpolation from the array.
//========================================================================================
{
	int   	LI;		// 	Lower index
	double	dx;
	
	dx	= (x_max-x_min)/nx;
	LI	= (int)floor((x_star-x_min)/dx);
	
	//	safety
	if (LI <= 0 ) LI = 0;
	if (LI >= nx-1 ) LI = nx-1;

	return	k_intp_linear(x_min + dx*LI, x_min + dx*(LI+1), fx[LI], fx[LI+1], x_star);
}
