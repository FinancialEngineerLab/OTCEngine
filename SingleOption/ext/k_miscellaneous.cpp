
#include <k_miscellaneous.hpp>
#include <stdexcept>
#include <vector> 

double *** k_misc_new_double_array_HS( int I, int J, int K )
//	dynamic memory allocation 
//	for 3-dim double array
{
	int		i, j;
	
	double	***return_address;
	
	return_address = new double** [I];
	
	for ( i = 0; i < I; i++ ) {
		return_address[i] = new double* [J];
		for ( j = 0; j < J; j++ )
			return_address[i][j] = new double [K];
	}
	return return_address;
}

double ** k_misc_new_double_array_HS( int I, int J )
//	dynamic memory allocation 
//	for 2-dim double array
{

	int		i;

	double	**return_address;

	return_address = new double* [I];

	for ( i = 0; i < I; i++ )
		return_address[i] = new double [J];

	return return_address;
}

int ** k_misc_new_int_array_HS( int I, int J )
//	dynamic memory allocation 
//	for 2-dim double array
{
	int		i;
	
	int		**return_address;
	
	return_address = new int* [I];
	
	for ( i = 0; i < I; i++ )
		return_address[i] = new int [J];

	return return_address;
}

void 	k_misc_copy_array_HS( int *target, int *source, int I )
{
	int		i;
	
	for ( i = 0; i < I; i++ )
		target[i] = source[i];
}

void 	k_misc_copy_array_HS( double *target, double *source, int I )
{
	int		i;
	
	for ( i = 0; i < I; i++ )
		target[i] = source[i];
}

void 	k_misc_copy_array_HS( double **target, double **source, int I, int J )
{
	int		i, j;
	
	for ( i = 0; i < I; i++ ) for ( j = 0; j < J; j++ )
		target[i][j] = source[i][j];
}

void 	k_misc_swap_address_HS( double ****A, double ****B )
//	swap address of double** *(e.g. 3-dim double array)
{
	double ***C;
	C = *A;
	*A = *B;
	*B = C;
}

void 	k_misc_swap_address_HS( double ***A, double ***B )
//	swap address of double** (e.g. 2-dim double array)
{
	double **C;
	C = *A;
	*A = *B;
	*B = C;
}

void 	k_misc_swap_address_HS( double **A, double **B )
//	swap address of double*	(e.g. 1-dim double array)
{
	double *C;
	C = *A;
	*A = *B;
	*B = C;
}

void	k_misc_delete_array_HS( double ***A, int I, int J )
//	clear dynamic memory allocation 
//	for 3-dim double array: 	A[I][J][K]
{
	int 	i, j;
	
	if ( A ) {
		for ( i = 0; i < I; i++ ) {
			if ( A[i] ) {
				for ( j = 0; j < J; j++ )
					if ( A[i][j] ) delete[] A[i][j];
				delete A[i];	
			}	
		}
		delete[] A;
	}
	A = NULL;
}

void	k_misc_delete_array_HS( double **A, int I )
//	clear dynamic memory allocation 
//	for 2-dim double array:	A[I][J]
{
	int 	i;
	
	if ( A ) {
		for ( i = 0; i < I; i++ )
			if ( A[i] ) delete[] A[i];
		delete[] A;
	}
	A = NULL;
}

void	k_misc_delete_array_HS( int **A, int I )
//	clear dynamic memory allocation 
//	for 2-dim int array:	A[I][J]
{
	int 	i;
	
	if ( A ) {
		for ( i = 0; i < I; i++ )
			if ( A[i] ) delete[] A[i];
		delete[] A;
	}
	A = NULL;
}

void	k_misc_delete_array_HS( double *A )
//	clear dynamic memory allocation 
//	for 1-dim double array:	A[I]
{
	if ( A )
		delete[] A;
	A = NULL;
}

void	k_misc_delete_array_HS( int *A )
//	clear dynamic memory allocation 
//	for 1-dim int array:	A[I]
{
	if ( A )
		delete[] A;
	A = NULL;
}

int k_misc_binsearch_lower_index( double *sorted_array, int num_ele, double key )
//	sorted_array is of increasing order
//	returns LI s.t. sorted_array[LI]<= key < sorted_array[LI+1]
//	return -1		if	key < sorted_array[0]
//	return num_ele	if	sorted_array[num_ele-1] <= key 
//	return -2		if  search failed
{

	if ( key < sorted_array[0] )
		return -1;
	if ( sorted_array[num_ele-1] <= key )
		return num_ele;

	int first=0,last=num_ele-1;
	int mid;

	while(first<=last)
	{
		mid=(first+last)/2;
		if( sorted_array[mid] <= key && key < sorted_array[mid+1] )
			return mid;
		if(key<sorted_array[mid])
			last=mid-1;
		else
			first=mid+1;
	}             
	return -2;
}

int check_strict_increasing( int *a, int nb_a )
//	return 0 if the array is strictly increasing
//	return 1 otherwise
{
	int	i;
	for(i=0; i<nb_a-1; i++){
		if(a[i] >= a[i+1]) return 1;
	}
	return 0;
}

int check_monotone_increasing( int *a, int nb_a )
//	return 0 if the array is monotone increasing
//	return 1 otherwise
{
	int	i;
	for(i=0; i<nb_a-1; i++){
		if(a[i] > a[i+1]) return 1;
	}
	return 0;
}


void xlvolToIvol(double** Ivol,double* Ivolsurface_xl, long rows, long cols, long startindex)
//skip first startindex of Ivol
//e.g, Ivolsurface(0 to M, 0 to N) with startindex=1, Ivol takes IvolSurface(1 to rows, 1 to cols)
{
    for(long i=0;i<rows;i++)
    {
        for(long j=0;j<cols;j++)
        {
            Ivol[i][j]=*(Ivolsurface_xl+(i+startindex)+(j+startindex)*(rows+1));
        }
    }

}
void ivlv(double spot, double **Lvol_kospi,
               double **Ivol_kospi,
               double* vol_term, int nb_vol_term, double* vol_strike, int nb_vol_strike,
               double* r_curve, double* r_term,int nb_r,
               double* div_curve, double* div_term, int nb_div)
{
//    myfunction fct;
    double** dvol_dt;
    double** dvol_dstrike;
    double** dvol2_dstrike2;
    dvol_dt=new double*[nb_vol_term];
    dvol_dstrike=new double*[nb_vol_term];
    dvol2_dstrike2=new double*[nb_vol_term];

    for(int i=0;i<nb_vol_term;i++){
        dvol_dt[i]=new double[nb_vol_strike];
        dvol_dstrike[i]=new double[nb_vol_strike];
        dvol2_dstrike2[i]=new double[nb_vol_strike];
    }


//    double dvol_dt[10][17];
//    double dvol_dstrike[10][17];
//    double dvol2_dstrike2[10][17];
    //dv/dt
    for(int i=0;i<nb_vol_term;i++){
        if(i==0){
            for(int j=0;j<nb_vol_strike;j++)
                dvol_dt[i][j]=(Ivol_kospi[i+1][j]-Ivol_kospi[i][j])/(vol_term[i+1]-vol_term[i]);
        }else if(i==(nb_vol_term-1)){
            for(int j=0;j<nb_vol_strike;j++)
                dvol_dt[i][j]=(Ivol_kospi[i][j]-Ivol_kospi[i-1][j])/(vol_term[i]-vol_term[i-1]);
        }else{
            for(int j=0;j<nb_vol_strike;j++)
                dvol_dt[i][j]=(Ivol_kospi[i+1][j]-Ivol_kospi[i-1][j])/(vol_term[i+1]-vol_term[i-1]);
        }
    }

    //dv/dk
    for(int j=0;j<nb_vol_strike;j++){
        if(j==0){
            for(int i=0;i<nb_vol_term;i++){
                dvol_dstrike[i][j]=(Ivol_kospi[i][j+1]-Ivol_kospi[i][j])/(vol_strike[j+1]-vol_strike[j]);
                dvol2_dstrike2[i][j]=(Ivol_kospi[i][j+2]-2.0*Ivol_kospi[i][j+1]+Ivol_kospi[i][j])/(vol_strike[j+2]-vol_strike[j+1])/(vol_strike[j+1]-vol_strike[j]);
            }
        }else if(j==nb_vol_strike-1){
            for(int i=0;i<nb_vol_term;i++){
                dvol_dstrike[i][j]=(Ivol_kospi[i][j]-Ivol_kospi[i][j-1])/(vol_strike[j]-vol_strike[j-1]);
                dvol2_dstrike2[i][j]=(Ivol_kospi[i][j]-2.0*Ivol_kospi[i][j-1]+Ivol_kospi[i][j-2])/(vol_strike[j]-vol_strike[j-1])/(vol_strike[j-1]-vol_strike[j-2]);
            }
        }else{
            for(int i=0;i<nb_vol_term;i++){
                dvol_dstrike[i][j]=(Ivol_kospi[i][j+1]-Ivol_kospi[i][j-1])/(vol_strike[j+1]-vol_strike[j-1]);
                dvol2_dstrike2[i][j]=(Ivol_kospi[i][j+1]-2.0*Ivol_kospi[i][j]+Ivol_kospi[i][j-1])/(vol_strike[j+1]-vol_strike[j])/(vol_strike[j]-vol_strike[j-1]);
            }
        }
    }

    for(int i=0;i<nb_vol_term;i++)
        for(int j=0;j<nb_vol_strike;j++){
            double v, dvdk,dvdt,dv2dk2;
            double d1,d2,t,k,r,q;
            v=Ivol_kospi[i][j];
            dvdk=dvol_dstrike[i][j];
            dvdt=dvol_dt[i][j];
            dv2dk2=dvol2_dstrike2[i][j];
            t=vol_term[i];
            k=vol_strike[j];
            r=intp1d(r_term,r_curve,nb_r,t);
            q=intp1d(div_term,div_curve,nb_div,t);
            d1=(log(spot/k)+(r-q+0.5*v*v)*t)/(v*sqrt(t));
            d2=d1-v*sqrt(t);
            Lvol_kospi[i][j]=sqrt(v*v+2.0*t*v*dvdt+2*(r-q)*k*t*v*dvdk)/sqrt(1.0+2.0*k*d1*sqrt(t)*dvdk+k*k*t*(d1*d2*dvdk*dvdk+v*dv2dk2));

        }

    for(int i=0;i<nb_vol_term;i++){
        delete[] dvol_dt[i];
        delete[] dvol_dstrike[i];
        delete[] dvol2_dstrike2[i];
    }

    delete[] dvol_dt;
    delete[] dvol_dstrike;
    delete[] dvol2_dstrike2;
}


double intp1d(double* tenor, double* curve, int nCol, double t)
{
	double result;

	if (t<=tenor[0])
	{
		t=tenor[0];
		//result=curve[0];
		//return result;
	}
    else if (t>=tenor[(long)nCol-1])
	{
		t=tenor[(long)nCol-1];
		//result=curve[(long)nCol-1];
		//return result;
	}

    for (long i=1; i<(long)nCol; i++)
	{
        if (t <= tenor[i])
		{
            result = (curve[i] * (t - tenor[i - 1]) + curve[i - 1] * (tenor[i] - t)) / (tenor[i] - tenor[i - 1]);
            return result;
		}
	}
	throw std::logic_error("interpolation fails");
	return 1;
}

//1d interpolation 
double intp1d(double target, double *x, double *y, int min, int max)
{
	if(target <= x[min])
		target=x[min];
	else if(target > x[max])
		target=x[max];
	for(int i=min+1;i<=max;i++){
		if(target <= x[i])
			return (y[i]*(target-x[i-1])+y[i-1]*(x[i]-target))/(x[i]-x[i-1]);
	}

	throw std::logic_error("interpolation fails");
	return 1;
}

double intp1d(double target, double* px, std::vector<double>& uold, int min, int max, int init_i)
{
	if (target <= px[min])
		target = px[min];
	else if (target > px[max])
		target = px[max];
	for (int i = min + 1; i <= max; i++) {
		if (target <= px[i])
			return (uold[i] * (target - px[i - 1]) + uold[i - 1] * (px[i] - target)) / (px[i] - px[i - 1]);
	}

	throw std::logic_error("interpolation fails");
	return 1;
}
double intp1d_delta(double target, double *x, double *y, int min, int max)
{
	if (target <= x[min + 1]) 
		target = x[min + 1];
	else if (target >= x[max-1])
		target = x[max - 1];

	for (int i = min + 1; i <= max; i++) {
		if (target <= x[i]) {
			return (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
		}
	}

	throw std::logic_error("interpolation fails");
	return 1;
}

double getforward(double* tenor, double* curve, int nCol, double t)
{
	double result, drdt;

	if (t<=tenor[0])
	{
		result=curve[0];
		return result;
	}
    else if (t>=tenor[(long)nCol-1])
	{
		result=curve[(long)nCol-1];
		return result;
	}

    for (long i=0; i<(long)nCol-1; i++)
	{
        if (t < tenor[i + 1] && t >= tenor[i])
		{
            drdt = (curve[i + 1] - curve[i]) / (tenor[i + 1] - tenor[i]);
            result = ((t - tenor[i]) * curve[i + 1] + (tenor[i + 1] - t) * curve[i]) / (tenor[i + 1] - tenor[i]);
            result = result + t * drdt;
			return result;
		}
	}

	throw std::logic_error("interpolation fails in getforward()");
	return 1;

}

double getforward(const double tau, const double *r, const double *term, const long nb)
{
	if(tau<=term[0]) return r[0];
	if(tau>=term[nb-1]) return r[nb-1];
	for(long i=0;i<nb;i++)
	{
		if(tau <term[i+1] && tau >=term[i]){
			double drdt;
			drdt=(r[i+1]-r[i])/(term[i+1]-term[i]);
			return tau*drdt+((tau - term[i]) * r[i+1] + (term[i+1] - tau) * r[i]) / (term[i+1] - term[i]);  //forward rate=r(t)+dt/dt*t
		}
	}
	throw std::logic_error("interpolaton fail :getforward");
	return 1;
}

string getFnameTimeStartingWith(string init_str)
{
	std::ostringstream oss;


	time_t curr_time;
	struct tm *curr_tm;

	curr_time = time(NULL);
	curr_tm = localtime(&curr_time);
	string str_mon;
	string str_day;
	string str_hour;
	string str_min;

	if (curr_tm->tm_mon + 1 < 10)
		str_mon = string("0") + to_string(curr_tm->tm_mon + 1);
	else
		str_mon = to_string(curr_tm->tm_mon + 1);

	if (curr_tm->tm_mday< 10)
		str_day = string("0") + to_string(curr_tm->tm_mday);
	else
		str_day = to_string(curr_tm->tm_mday);

	if (curr_tm->tm_hour< 10)
		str_hour = string("0") + to_string(curr_tm->tm_hour);
	else
		str_hour = to_string(curr_tm->tm_hour);

	if (curr_tm->tm_min  < 10)
		str_min = string("0") + to_string(curr_tm->tm_min);
	else
		str_min = to_string(curr_tm->tm_min);


	oss << init_str << str_mon << str_day << str_hour << str_min << ".csv";

	return oss.str();
}
