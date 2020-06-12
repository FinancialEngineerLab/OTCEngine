#ifndef KDB_LINEARALGEBRA_H
#define KDB_LINEARALGEBRA_H

#include <vector>
#include <cmath>
#include <limits>
#include <iostream> //cout
using namespace std;
/**********************************************************************/
/* Program Name      : Linear Algebra library
/* version           : 0.0.0.0
/* author            : kicheon chang
/* date              : 2008.12.26
/* Modified by       : kicheon chang
/* Modified at       : 
/* Copyright         : KDB Quant team
/* Description       : Matrix operation, Linear system solver, etc
/* Related Doc. Name : none
/**********************************************************************/

typedef double DP;

/************************************************************************/
// Class for matrix
//
/************************************************************************/
template <class Object> class matrix // Object is just a name by programmer
{
public: 
	matrix(int rows=0, int cols=0):array(rows) // you can give a size in advance or later
	{
		for (int i=0;i<rows;i++)
			array[i].resize(cols);
	}
	int numrows() const		{ return (int)array.size(); }
	int numcols() const		{ return numrows() >0 ? (int)array[0].size(): 0;}
	void resize(int rows,int cols){
		array.resize(rows);
		for(int i=0;i<rows;i++)
			array[i].resize(cols);
	}
	void operator=(const matrix<double> &a);
	const vector<Object> & operator[] (int row) const 		{ return array[row];}
	vector<Object> & operator[] (int row)		{ return array[row];}

private:
	vector < vector<Object> > array;
};

//FUNCTION PROTOTYPE
//operator overloading...
matrix<double> operator*(const matrix<double> &a, const matrix<double> &b);
matrix<double> operator*(const matrix<double> &a, double b);
matrix<double> operator*(double b, const matrix<double> &a);
vector<double> operator*(const matrix<double> &a, const vector<double> &b);
matrix<double> operator/(const matrix<double> &a, double b);
matrix<double> operator+(const matrix<double> &a, const matrix<double> &b);
matrix<double> operator-(const matrix<double> &a, const matrix<double> &b);
// "= " is already overloaded in vector class

matrix<double> k_mx_transpose(matrix<double> & m);
matrix<double> k_mx_inverse(matrix<double> & mx);
matrix<double> k_mx_lu_decomp(matrix<double> & mx);
matrix<double> k_mx_chole_decomp(matrix<double> & mx);
void Chol_Decomp(double **A, double **L, int n, double eps);


matrix<double> k_mx_unit_matrix(int row_size, int col_size);
matrix<double> k_mx_null_matrix(int row_size, int col_size);
matrix<double> k_mx_identity_matrix(int row_col_size);

int	k_mx_tridiagonal_solver(int array_size, double *a,double* b, double* c, double *y, double *x);
int	k_mx_tridiagonal_solver(int array_size, double *coef, double *y, double *x);
int k_mx_tridiagonal_solver(double* a,double* b,double* c, double* y, int start_index_num,int end_index_num);
int k_mx_tridiagonal_solver(vector<double> a,vector<double>b,vector<double>c, vector<double> y, vector<double> &x);
int k_mx_linear_system_solver(matrix<double> A,vector<double> b, vector<double> &x);

//	added by HANGSEOB
double k_mx_Alpha_System_Solver(int dim_a,
							  double **alpha,    //    a square-matriix alpha
							  double *x,                //   output :  solution vector 
							  double *beta              //   vector                          [alpha]x = beta
							  );


void k_mx_cout(matrix<double> &mx);
void k_mx_cout(vector<double> &vt);
double k_mx_determinant(matrix<double> &mx);
matrix<double> k_mx_row_operation(matrix<double> &mx,double a,int row_num);
matrix<double> k_mx_col_operation(matrix<double> &mx,double a,int col_num);
matrix<double> k_mx_row_exchange(matrix<double> &mx, int row_num01, int row_num02);
matrix<double> k_mx_col_exchange(matrix<double> &mx, int col_num01, int col_num02);
vector<double> k_mx_column_extract(matrix<double> &mx, int col_num);

void k_mx_jacobi_eigenvalue(matrix<double> &a, int iter,matrix<double> &returnmat);
void k_mx_jacobi_eigenvector(matrix<double> &mat, int iter, matrix<double> &returnmat);

//Internal use only
inline void		rot(matrix<double> &a, const DP s, const DP tau, const int i,const int j, const int k, const int l);
static double	norm_inf(matrix<double> &x);
static double	norm2(matrix<double> &x);
static int		find_p(matrix<double> &x);
static void		find_max(matrix<double> &x, int &p, int &q);

void mx_powermethod(matrix<double> a, vector<double> &out_vector, double &out_value,double epsilon, int iter);
//this is not varified!! do not use!!
#endif