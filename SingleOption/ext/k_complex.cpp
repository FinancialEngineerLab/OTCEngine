#include <kdb_complex.hpp>


double k_com_abs(k_complex z)  /* calculate the absolute value of z */
{ 
   return(sqrt(z.re*z.re + z.im*z.im));
}

k_complex k_com_rmul(double r, k_complex z)
{
   k_complex val;

   val.re = r*z.re;
   val.im = r*z.im;

   return(val);
}


k_complex k_com_exp(k_complex z)
{
	k_complex val;

	val.re = exp(z.re)*cos(z.im);
	val.im = exp(z.re)*sin(z.im);

	return(val);
}



k_complex k_com_mul(k_complex z1,k_complex z2)  /* calculate z1 * z2 */
{
  k_complex z3; 
  z3.re = z1.re * z2.re - z1.im * z2.im;
  z3.im = z1.re * z2.im + z1.im * z2.re;
  return(z3);
}

k_complex k_com_add(k_complex z1,k_complex z2)  /* calculate z1 + z2 */
{
   k_complex z3; 
  z3.re = z1.re + z2.re;
  z3.im = z1.im + z2.im;

  return(z3);
}

k_complex k_com_sub(k_complex z1,k_complex z2)  /* calculate z1 - z2 */
{
  k_complex z3; 
  z3.re = z1.re - z2.re;
  z3.im = z1.im - z2.im;
  return(z3);
}

k_complex k_com_div(k_complex z1,k_complex z2)  /* calculate z1 / z2 */
{ 
  k_complex z3;
  z3.re = (z1.re * z2.re + z1.im * z2.im)/(k_com_abs(z2)*k_com_abs(z2));
  z3.im = (-z1.re * z2.im + z1.im * z2.re)/(k_com_abs(z2)*k_com_abs(z2));
  return(z3);
}


k_complex k_com_sqrt(k_complex z) /* calculate sqrt(z) */
{ 
   k_complex z1;
   double theta,r;
 
   theta = 0.5*atan2(z.im,z.re);
   r = sqrt(k_com_abs(z));
   z1.re = r*cos(theta); z1.im = r*sin(theta); 
  
   return(z1); 
} 

k_complex k_com_inv(k_complex z)
{
	k_complex id;
	
	id.re = 1.0;
	id.im = 0.0;

    return(k_com_div(id,z));
}


k_complex k_com_log(k_complex z)
{
	k_complex tp;

	double r;
	double theta;

	r = sqrt(z.re*z.re + z.im*z.im);

   theta = atan2(z.im, z.re);

   tp.re = log(r);
   tp.im = theta;

   return(tp);
}