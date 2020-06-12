#include <kdb_LinearAlgebra.hpp>

//FUNCTION DEFINITION...

/************************************************************************/
// Operator(=)
//
/************************************************************************/
void matrix<double>::operator =(const matrix<double> &a){
	int p =a.numrows(), q=a.numcols();
	this->resize(p,q);
	for(int i=0;i<p;i++) for(int j=0;j<q;j++){
		(*this)[i][j]=a[i][j];
	}
}

/************************************************************************/
// Operator(*)
// matrix*matrix
/************************************************************************/
matrix<double> operator*(const matrix<double> &a, const matrix<double> &b)
{
	int p =a.numrows(), q=a.numcols(), r=b.numrows(), s=b.numcols();
	matrix<double> c(p,s); // returned matix could not be square-matrix... 
	if (q!=r) //cannot define multiply
	{
		c[0][0]=0.0;
	} //do nothing
	else
	{
		//initialized returned matrix
		for (int i=0;i<p;i++)	
			for (int j=0;j<s;j++)
				c[i][j]=0;
		//calc multiply.....
		for (int i=0;i<p;i++) //row index
			for (int j=0 ; j<s;j++) //column index, check that it returns p*s matrix
				for (int k=0 ; k<q ;k++) //calc elements
					c[i][j] += a[i][k] * b[k][j];
	}
	return c;
}

/************************************************************************/
// Operator(+)
// matrix+matrix
/************************************************************************/
matrix<double> operator+(const matrix<double> &a, const matrix<double> &b)
{
	int p =a.numrows(), q=a.numcols(), r=b.numrows(), s=b.numcols();
	matrix<double> c(p,s); // returned matix could not be square-matrix... 
	if (p==r && q==s) //cannot define multiply
	{
		for (int i=0;i<p;i++)
			for (int j=0;j<s;j++)
				c[i][j]=0;
		for (int i=0;i<p;i++)
			for (int j=0 ; j<s;j++)
				c[i][j] += a[i][j] + b[i][j];
		return c;
	}
	else
	{
		c[0][0]=0.0;
	} //do nothing
	return c;
}

/************************************************************************/
// Operator(-)
// matrix-matrix
/************************************************************************/
matrix<double> operator-(const matrix<double> &a, const matrix<double> &b)
{
	int p =a.numrows(), q=a.numcols(), r=b.numrows(), s=b.numcols();
	matrix<double> c(p,s); // returned matix could not be square-matrix... 
	if (p==r && q==s) //cannot define multiply
	{
		for (int i=0;i<p;i++)
			for (int j=0;j<s;j++)
				c[i][j]=0;
		for (int i=0;i<p;i++)
			for (int j=0 ; j<s;j++)
				c[i][j] += a[i][j] - b[i][j];
		return c;
	}
	else
	{
		c[0][0]=0.0;
	} //do nothing
	return c;
}
/************************************************************************/
// Operator(*)
// matrix*scalar
/************************************************************************/
matrix<double> operator*(const matrix<double> &a, double b)
{
	int p =a.numrows(), q=a.numcols();
	matrix<double> c(p,q); 

	for (int i=0;i<p;i++)
		for (int j=0;j<q;j++)
			c[i][j]=0;
	for (int i=0;i<p;i++)
		for (int j=0 ; j<q;j++)
			c[i][j] = a[i][j]*b;
	return c;
}

/************************************************************************/
// Operator(*)
// scalar*matrix
/************************************************************************/
matrix<double> operator*(double b, const matrix<double> &a)
{
	int p =a.numrows(), q=a.numcols();
	matrix<double> c(p,q); 

	for (int i=0;i<p;i++)
		for (int j=0;j<q;j++)
			c[i][j]=0;
	for (int i=0;i<p;i++)
		for (int j=0 ; j<q;j++)
			c[i][j] = a[i][j]*b;
	return c;
}
/************************************************************************/
// Operator(*)
// matrix*vector; caution! vector*matrix is not defined.
/************************************************************************/
vector<double> operator*(const matrix<double> &a, const vector<double> &b)
{
	int i, n;
	n=(int)b.size();
	matrix<double> tmp(n,1);
	for(i=0;i<n;i++)
		tmp[i][0]=b[i];

	vector<double> tmp_vector(n);
	tmp=a*tmp;
	for(i=0;i<n;i++)
		tmp_vector[i]=tmp[i][0];
	return tmp_vector;
}

/************************************************************************/
// Operator(/)
// matrix/scalar
/************************************************************************/
matrix<double> operator/(const matrix<double> &a, double b)
{
	int p =a.numrows(), q=a.numcols();
	matrix<double> c(p,q); 

	for (int i=0;i<p;i++)
		for (int j=0;j<q;j++)
			c[i][j]=0;
	for (int i=0;i<p;i++)
		for (int j=0 ; j<q;j++)
			c[i][j] = a[i][j]/b;
	return c;
}
/************************************************************************/
// cout utility for matrix
//
/************************************************************************/
void k_mx_cout(matrix<double> &mx)
{
	int n=mx.numrows();
	int m=mx.numcols();
	int i,j;
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			cout<< mx[i][j] << '\t' ;
		}
		cout << endl;
	}
}
/************************************************************************/
// cout utility for vector
//
/************************************************************************/
void k_mx_cout(vector<double> &vt)
{
	int n=(int)vt.size();
	for(int i=0;i<n;i++)
		cout << vt[i] << endl;
}
/************************************************************************/
// Determinant
//
/************************************************************************/
double k_mx_determinant(matrix<double> &mx)
{
	int i;
	double tmp_prd=1.0;
	matrix<double> a;
	a=k_mx_lu_decomp(mx);
	for(i=0;i<a.numcols();i++)
	{
		tmp_prd*=a[i][i];
	}
	return tmp_prd;
}
/************************************************************************/
// LU decomposition
//
/************************************************************************/
matrix<double> k_mx_lu_decomp(matrix<double> &mx)
{
	int i, j,k;
	double tmpSum;
	matrix<double> a(mx.numrows(),mx.numcols());
	matrix<double> b(mx.numrows(),mx.numcols());
	matrix<double> c(mx.numrows(),mx.numcols());

	for(j=0;j<mx.numcols();j++)
	{
		for(i=0;i<=j;i++)
		{
			tmpSum=0.0;
			for(k=0;k<=i-1;k++)
			{
				tmpSum+=a[i][k]*b[k][j];
			}
			if(j==0) tmpSum=0.0;
			b[i][j]=mx[i][j]-tmpSum;
			c[i][j]=b[i][j];
		}
		for(i=j+1;i<mx.numrows();i++)
		{
			tmpSum=0.0;
			for(k=0;k<=j-1;k++)
			{
				tmpSum+=a[i][k]*b[k][j];
			}
			if(i==0) tmpSum=0.0;
			a[i][j]=1.0/b[j][j]*(mx[i][j]-tmpSum);
			c[i][j]=a[i][j];
		}
	}
	return c;
}
/************************************************************************/
// tridiangonal equation solver
//
/************************************************************************/
int k_mx_tridiagonal_solver(vector<double> a,vector<double>b,vector<double>c, vector<double> y, vector<double> &x)
{
	int i;
	int n=(int)a.size();
	//we accept n-elements vectors; a, b, c, y, x
	//do not use a[0] and c[n-1]
	//remove a[i]'s
	for(i=1;i<n;i++){
		b[i]=b[i]-c[i-1]*a[i]/b[i-1];
		y[i]=y[i]-y[i-1]*a[i]/b[i-1];
	}
	//remove c[i]'s
	for(i=n-2;i>=0;i--)
		y[i]=y[i]-y[i+1]*c[i]/b[i+1];

	//return solutions
	for(i=0;i<n;i++)
		x[i]=y[i]/b[i];

	return 1;
}

/************************************************************************/
// tridiangonal equation solver using array instead of vector
// the original interface a,b,c,y,x are the vectors with array_size elements
/************************************************************************/
int	k_mx_tridiagonal_solver(int array_size, double *a,double* b, double* c, double *y, double *x)
{
	int i,n=array_size;
	//remove a[i]'s
	for(i=1;i<n;i++){
		b[i]=b[i]-c[i-1]*a[i]/b[i-1];
		y[i]=y[i]-y[i-1]*a[i]/b[i-1];
	}
	//remove c[i]'s
	for(i=n-2;i>=0;i--)
		y[i]=y[i]-y[i+1]*c[i]/b[i+1];

	//return solutions
	for(i=0;i<n;i++)
		x[i]=y[i]/b[i];

	return 1;
}
/************************************************************************/
// tridiangonal equation solver using array instead of vector
// a,b,c are coefficients
// array index can be selected by user; y is b vector of Ax=b, also acts as returnig container 
//
//=================== Hangseob added more comments at 090603 ======================
//	tri-diagonal 행렬 방정식을 푼다.
//	방정식	a[i]V[i-1] + b[i]V[i] + c[i]V[i+1] = y[i], i=M,2,3,...,N을 만족하는 
//	V[M], V[M+1], ... , V[N]을 구한다. 
//		(N=end_index_num,M=start_index_num) 
//
//	a[M] = c[N] = 0을 만족할 때이고 결과는 y[M], y[M+1], ... , y[N]에 저장된다.
//	(a[M] = c[N] = 0이 아니더라도 a[M] = c[N] = 0이라 가정하고 푼다. )
//
//	a, b, c 은 값만 받고처리하고 original value는 그대도 보존한다. 
//
//
/************************************************************************/
int k_mx_tridiagonal_solver(double* a,double* b,double* c, double* y, int start_index_num,int end_index_num)
{
	int N=end_index_num;
	int M=start_index_num;
	int i;
	double *A, *C, *B;	

	A = (double *)malloc((N+1)*sizeof(double));
	B = (double *)malloc((N+1)*sizeof(double));
	C = (double *)malloc((N+1)*sizeof(double));

	for (i=M;i<=N;i++) {
		A[i] = a[i];
		B[i] = b[i];
		C[i] = c[i];
	}

	A[M] = C[N] = 0.0;

	if (M==N) y[M] /= B[M];
	else {
		for (i=M+1;i<=N;i++) {
			A[i] /= B[i-1];
			B[i] -= A[i]*C[i-1];
			y[i] -= A[i]*y[i-1];
		}
		y[N] /= B[N];

		for (i=N-1;i>=M;i--) y[i] = (y[i]-C[i]*y[i+1])/B[i];
	}

	if (A) free(A);
	if (B) free(B);
	if (C) free(C);

	return 1;
}

int	k_mx_tridiagonal_solver(int array_size, double *coef, double *y, double *x)
{
	int i,n=array_size;
	double *a;
	double *b;
	double *c;

	a=(double*)calloc(n+1,sizeof(double));
	b=(double*)calloc(n,sizeof(double));
	c=(double*)calloc(n+1,sizeof(double));

	a[0]=0.;
	b[0]=coef[0];
	c[0]=coef[1];
	for(i=1;i<n-1;i++){
		a[i]=coef[3*i-1];
		b[i]=coef[3*i+0];
		c[i]=coef[3*i+1];}
	a[n-1]=coef[3*n-4];
	b[n-1]=coef[3*n-3];
	c[n-1]=0.;

	i=k_mx_tridiagonal_solver(n,a,b,c,y,x);

	if(a) 		free(a);
	if(b) 		free(b);
	if(c) 		free(c);

	return 1;
}


/************************************************************************/
// matrix transpose
//
/************************************************************************/
matrix<double> k_mx_transpose(matrix<double> & m)
{
	int p=m.numrows(), q=m.numcols();
	int i,j;
	matrix<double> c(q,p);
	for(i=0;i<q;i++)
	{	for(j=0;j<p;j++)
	{c[i][j]=m[j][i];}
	}
	return c;
}

/************************************************************************/
// finding inverse matrix
//
/************************************************************************/
matrix<double> k_mx_inverse(matrix<double> & mx)
{
	int p=mx.numrows(), q=mx.numcols();
	int i,j,k;
	//matrix<double> mx(q,p);
	matrix<double> oldmx(q,p);
	matrix<double> oldix(q,p);
	matrix<double> ix(q,p);
	//mx=m;
	//double a,b,c,d,e,f,g,h,tmp;		
	if(p!=q) return ix;

	//make identity matrix
	for(i=0;i<q;i++){
		for(j=0;j<q;j++){
			if(j==i) ix[i][j]=1.0;
			else ix[i][j]=0.0;
		}
	}

	//set non-diagonal elements zero
	for(j=0;j<q;j++){//for j-th column...
		for(i=0;i<q;i++){ //for i-th elements
			if(i!=j){
				oldmx=mx;
				oldix=ix;
				for(k=0;k<q;k++){
					mx[i][k]=oldmx[i][k]-oldmx[j][k]/oldmx[j][j]*oldmx[i][j];
					ix[i][k]=oldix[i][k]-oldix[j][k]/oldmx[j][j]*oldmx[i][j];
				}
			}
		}
	}

	//make all the diangoanl elements "1"
	for(i=0;i<q;i++){
		for(j=0;j<q;j++){
			ix[i][j]=ix[i][j]/mx[i][i];
		}
	}
	return ix;
}
/************************************************************************/
// finding cholesky-decomposed matrix
//
/************************************************************************/
matrix<double> k_mx_chole_decomp(matrix<double> & mat)
{
	int N=mat.numcols(); //[row][column]
	int i,j,k;
	double s;
	matrix<double> L(N,N);
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			L[i][j]=0.0;
		}
	}
	for(j=0;j<N;j++){
		s=0.0;
		for(k=0;k<j;k++){
			s+=L[j][k]*L[j][k];
		}
		L[j][j]=mat[j][j]-s;
		L[j][j]=sqrt(L[j][j]);
		for(i=j+1;i<N;i++){
			s=0.0;
			for(k=0;k<j;k++){
				s+=L[i][k]*L[j][k];
			}
			L[i][j]=(mat[i][j]-s)/L[j][j];
		}
	}
	return L;
}

void Chol_Decomp(double **A, double **L, int n, double eps)
{
	int i,j,k;
	double sum;

	for (i=0;i<n;i++) 
	{
		for (j=0;j<n;j++) 
		{
			L[i][j] = A[i][j];
		}
	}

	for (i=0;i<n;i++) 
	{
		for (j=0;j<n;j++) 
		{
			sum = A[i][j];
			for (k=i-1;k>=0;k--) sum -= L[i][k]*L[j][k];
			if (i==j) 
			{
				if (sum<=0.0) sum = eps;
				L[i][i]=sqrt(sum);
			}
			else L[j][i]=sum/L[i][i];
		}
	}

	for (i=0;i<n;i++) 
	{
		for (j=i+1;j<n;j++) 
		{
			L[i][j] = 0.0;
		}
	}
}

/*************************************************************/
// functions for Eigen Value and Vectors
/*************************************************************/

/************************************************************************/
// finding cholesky-decomposed matrix
//
/************************************************************************/
inline void k_mx_rot(matrix<double> &a, const DP s, const DP tau, const int i,
					 const int j, const int k, const int l)
{
	DP g,h;

	g=a[i][j];
	h=a[k][l];
	a[i][j]=g-s*(h+g*tau);
	a[k][l]=h+s*(g-h*tau);
}

/************************************************************************/
// internal use only!!
//
/************************************************************************/
static double norm_inf(matrix<double> &x)
{
	int n=x.numcols(), i;
	double mx;
	mx=x[0][0];
	for(i=1;i<n;i++) mx=max(mx,x[0][i]);
	return mx;
}
/************************************************************************/
// internal use only!!
//
/************************************************************************/
static double norm2(matrix<double> &x)
{
	int n=x.numcols(), i;
	double sum=0.0;
	for(i=0;i<n;i++) sum+=x[0][i];
	return sum;
}
/************************************************************************/
// internal use only!!
//
/************************************************************************/
static int find_p(matrix<double> &x)
{
	int n=x.numcols(), i;
	for(i=0;i<n;i++) 
		if (fabs(x[0][i])==norm_inf(x)) return i;
	return 0;
}

/************************************************************************/
//find the dominant eigenvalue and an associated eigenvotor 
//this is not varified!! do not use!!
/************************************************************************/
void mx_powermethod(matrix<double> a, vector<double> &out_vector, double &out_value,double epsilon, int iter)
{
	int n=a.numcols(),p,i;

	double mu,error;
	matrix<double> x(1,n),y(1,n); //row vector
	for(i=0;i<n;i++)
		x[0][i]=1.0;

	p=find_p(x);
	x=x/double(p+1);
	for(i=0;i<iter;i++)
	{
		y=k_mx_transpose(a*k_mx_transpose(x));
		mu=y[0][p];
		if (mu==0.0) //error if yp==0
		{
			//eigenvector[0].resize(0);
		}
		p=find_p(y);
		error=norm_inf(x-(y/mu));
		x=y/y[0][p];
		if(error<epsilon) break;
	}
	out_vector=x[0];
	out_value=mu;
}
/************************************************************************/
// all eigenvalues of symmtric matrix using Jacobi rotation algo
//
/************************************************************************/
void k_mx_jacobi_eigenvalue(matrix<double> &mat, int iter, matrix<double> &returnmat)
{
	int p,q,i,j,k;
	int n=mat.numcols();
	double d,t,x,c,s; double check;
	matrix<double> u(n,n), w(n,n),b(n,n),a(n,n);
	a=mat;

	for(i=0;i<iter;i++)
	{
		for(j=0;j<n;j++) // set the size of return-matrix
			for(k=0;k<n;k++) 
			{
				u[j][k]=0.0; w[j][k]=0.0; b[j][k]=0.0;// a[i][j]=0.0;		
			}
			find_max(a, p, q);
			//check=a[q][q];
			//check=a[p][p];
			d=a[q][q]-a[p][p];
			if(d==0.0) {t=1.0;}
			else 
			{
				x=d/a[p][q]/2.0;
				t=x/fabs(x)  /  (fabs(x)+sqrt(x*x+1.0));
			}
			c=1.0/sqrt(t*t+1.0);
			s=t*c;
			for(j=0;j<n;j++)
			{
				u[j][j]=1.0;
				w[j][j]=1.0;
			}
			u[p][p]=c;	u[p][q]=s;
			u[q][p]=-s;	u[q][q]=c;
			w[p][p]=c;	w[p][q]=-s;
			w[q][p]=s;	w[q][q]=c;
			b=a*u;
			b=w*b;
			a=b;
			for(j=0;j<n;j++)
				for(k=0;k<n;k++)
				{check=a[j][k];}
				//cout << i<< endl ;
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i!=j) b[i][j]=0.0;
		}
	}
	returnmat=b;
}

/************************************************************************/
// find the maximum element in the matrix
//
/************************************************************************/
void find_max(matrix<double> &x, int &p, int &q)
{
	double big=0.0;
	int n=x.numcols();
	p=0; q=0;
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			if(i!=j && fabs(x[i][j])>big) 
			{
				big=fabs(x[i][j]);
				p=i,q=j;
			}
}
/************************************************************************/
//  approx eigenvectors of a symmetric matrix using Jacobi iterative algo
//
/************************************************************************/
void k_mx_jacobi_eigenvector(matrix<double> &mat, int iter, matrix<double> &returnmat)
{
	int n=mat.numcols();
	int i,j,k; //loop index
	int p,q;
	double d,t,x,c,s;
	matrix<double> a(n,n), v(n,n), b(n,n);
	matrix<double> u(n,n), w(n,n);
	a=mat;
	for(i=0;i<n;i++) // make identity matrix
		for(j=0;j<n;j++)
		{
			v[i][j]=0.0;
			b[i][j]=0.0;
			if(i==j) v[i][j]=1.0;
		}
		for(i=0;i<iter;i++)
		{
			find_max(a,p,q);
			d=a[q][q]-a[p][p];
			if(d==0.0) t=1.0;
			else
			{
				x=d/a[p][q]/2.0;
				t=x/fabs(x) /(fabs(x)+sqrt(x*x+1.0));
			}
			c=1.0/sqrt(t*t+1.0);
			s=t*c;
			for(j=0;j<n;j++)
				for(k=0;k<n;k++)
				{
					u[j][k]=0.0;
					w[j][k]=0.0;
				}
				for(j=0;j<n;j++)
				{
					u[j][j]=1.0;
					w[j][j]=1.0;
				}
				u[p][p]=c;	u[p][q]=s;
				u[q][p]=-s;	u[q][q]=c;
				w[p][p]=c;	w[p][q]=-s;
				w[q][p]=s;	w[q][q]=c;
				v=v*u;
				b=a*u;
				b=w*b;
				a=b;
		}
		returnmat=v;
}
/************************************************************************/
// generates m*n matrix of which elements are all unit
// m=row_size, n=col_size
/************************************************************************/
matrix<double> k_mx_unit_matrix(int row_size, int col_size)
{
	int i,j;
	matrix<double> return_mx(row_size,col_size);
	for(i=0;i<row_size;i++)
	{
		for(j=0;j<col_size;j++)
		{
			return_mx[i][j]=1.0;
		}
	}
	return return_mx;
}
/************************************************************************/
// generates m*n matrix of which elements are all null(zero)
// m=row_size, n=col_size
/************************************************************************/

matrix<double> k_mx_null_matrix(int row_size, int col_size)
{
	//int i,j;
	matrix<double> return_mx(row_size,col_size);
	return_mx=k_mx_unit_matrix(row_size,col_size)*0.0;
	return return_mx;
}
/************************************************************************/
// generates n*n identity matrix
// n=row_col_size
/************************************************************************/
matrix<double> k_mx_identity_matrix(int row_col_size)
{
	int i,j;
	matrix<double> return_mx(row_col_size,row_col_size);
	for(i=0;i<row_col_size;i++)
	{
		for(j=0;j<row_col_size;j++)
		{
			return_mx[i][j]=0.0;
			if(i==j) return_mx[i][j]=1.0;
		}
	}
	return return_mx;
}

/************************************************************************/
//	Linear system solver using a scaled pivoting	
//	find x  s.t  A*x = b; A is coefficient matrix, b is vector for value of function
/************************************************************************/
int k_mx_linear_system_solver(matrix<double> A,vector<double> b, vector<double> &x)
{
	int dim_a=(int)b.size();
	int i_true=1;
	int i_false=0;
	int 	i,j,k,mi;
	double temp,mult;
	vector<double> TA;
	vector<double> s(dim_a);
	x.resize(dim_a);

	for(j=0; j < dim_a; j++){
		s[j] = fabs(A[j][0]);
		for(k=1; k < dim_a; k++){
			if (s[j] < fabs(A[j][k]))
				s[j] = fabs(A[j][k]);
		}
		if (s[j] == 0.0){
			//no unique solution! //exit(0);
			return i_false;
		}
	}
	for(i=0; i < dim_a-1; i++){
		mi = i;
		for(k=i+1; k < dim_a; k++)
			if (fabs(A[k][i]/s[k]) > fabs(A[mi][i])/s[mi]){
				mi = k;
			}
			if (A[mi][i] == 0.0) {
				//pivot is zero! // exit(0);
				return i_false;
			}
			if (i != mi){
				TA = A[i];
				A[i] = A[mi];
				A[mi] = TA;

				temp = b[i];
				b[i] = b[mi];
				b[mi] = temp;
			}

			for(j=i+1; j < dim_a; j++) {
				mult = A[j][i] / A[i][i];
				for(k=i+1; k < dim_a; k++) {
					A[j][k] = A[j][k] - mult * A[i][k];
				}
				b[j] = b[j] - mult * b[i];
				A[j][i] = 0.0;
			}
	}

	// backward substitution procedure
	x[dim_a-1] = b[dim_a-1]/A[dim_a-1][dim_a-1];
	for(i = dim_a-2; i >= 0 ; i--) {
		temp = 0.0;
		for(j=i+1; j < dim_a; j++)
			temp = temp + (A[i][j] * x[j]);
		x[i] = (b[i] - temp) / A[i][i];
	}
	return i_true;
}

double k_mx_Alpha_System_Solver(
							  int dim_a,				//	dimension of alpha
							  double **alpha,			//  square-matrix;		use index from 0 to dim_a-1
							  double *x,              //  solution vector;	use index from 0 to dim_a-1
							  double *beta            //						use index from 0 to dim_a-1
							  )
//======================================================
//	solve a system of linear equation
//	
//	find x s.t
//	alpha*x = beta
//
//	use scaled pivotting
//
//	return 1.0 if it finds solution successfully
//=======================================================
{
	int 	i,j,k,mi;
	double temp,mult,*s;
	double **A, *b;
	double *TA;

	b = (double *)malloc(dim_a*sizeof(double));
	s = (double *)malloc(dim_a*sizeof(double));
	A = (double **)malloc(dim_a*sizeof(double *));
	for(i=0; i < dim_a; i++)
		A[i] = (double *)malloc(dim_a*sizeof(double));



	//   Keep alpha  & beta invariant.       Instead of [alpha]x = beta,  we use Ax =b 
	for(i=0; i < dim_a; i++){                   
		for(j=0; j < dim_a; j++){             
			A[i][j] = alpha[i][j];
		}
		b[i] = beta[i];
	}	                     

	//  Linear system solver using a scaled pivoting	   
	for(j=0; j < dim_a; j++){
		s[j] = fabs(A[j][0]);
		for(k=1; k < dim_a; k++){
			if (s[j] < fabs(A[j][k]))
				s[j] = fabs(A[j][k]);
		}
		if (s[j] == 0.0){
			return(0.0);
		}
	}
	for(i=0; i < dim_a-1; i++){
		mi = i;
		for(k=i+1; k < dim_a; k++)
			if (fabs(A[k][i]/s[k]) > fabs(A[mi][i])/s[mi]){
				mi = k;
			}
			if (A[mi][i] == 0.0) {
				return(0.0);
			}
			if (i != mi){
				TA = A[i];
				A[i] = A[mi];
				A[mi] = TA;

				temp = b[i];
				b[i] = b[mi];
				b[mi] = temp;
			}

			for(j=i+1; j < dim_a; j++) {
				mult = A[j][i] / A[i][i];
				for(k=i+1; k < dim_a; k++) {
					A[j][k] = A[j][k] - mult * A[i][k];
				}
				b[j] = b[j] - mult * b[i];
				A[j][i] = 0.0;
			}
	}

	// backward substitution procedure
	x[dim_a-1] = b[dim_a-1]/A[dim_a-1][dim_a-1];
	for(i = dim_a-2; i >= 0 ; i--) {
		temp = 0.0;
		for(j=i+1; j < dim_a; j++)
			temp = temp + (A[i][j] * x[j]);
		x[i] = (b[i] - temp) / A[i][i];
	}


	for(i=dim_a-1;i >= 0; i--)
		if (A[i]) free(A[i]);
	if (A) free(A);

	if (s) free(s);

	if (b) free(b);

	return(1.0);
}




/************************************************************************/
// returns matrix of which row_num-th row elements are multiplied by "a" 
// CAUTION!	row_num starts from "0" not "1"
/************************************************************************/
matrix<double> k_mx_row_operation(matrix<double> &mx,double a,int row_num)
{
	int m,n,i;
	n=mx.numrows();
	m=mx.numcols();
	matrix<double> tmp_mx(n,m);
	matrix<double> null_mx(n,m);
	tmp_mx=mx;
	null_mx=k_mx_null_matrix(n,m);
	if((row_num<0) || (row_num>=n))
	{
		return null_mx;
	}
	for(i=0;i<m;i++)
		tmp_mx[row_num][i]=a*mx[row_num][i];
	return tmp_mx;
}
/************************************************************************/
// returns matrix of which col_num-th column elements are multiplied by "a" 
// CAUTION!	row_num starts from "0" not "1"
/************************************************************************/
matrix<double> k_mx_col_operation(matrix<double> &mx,double a,int col_num)
{
	int m,n,i;
	n=mx.numrows();
	m=mx.numcols();
	matrix<double> tmp_mx(n,m);
	matrix<double> null_mx(n,m);
	tmp_mx=mx;
	null_mx=k_mx_null_matrix(n,m);
	if(col_num<0 || col_num>=m)	return null_mx;
	for(i=0;i<n;i++)
		tmp_mx[i][col_num]=a*mx[i][col_num];
	return tmp_mx;
}
/************************************************************************/
// returns matrix of which row_num01-th row is exchanged by row_num02-th row
// CAUTION!	row_num and col_num start from "0" not "1"
/************************************************************************/
matrix<double> k_mx_row_exchange(matrix<double> &mx, int row_num01, int row_num02)
{
	int m,n;//,i;
	n=mx.numrows();
	m=mx.numcols();
	matrix<double> tmp_mx(n,m);
	matrix<double> null_mx(n,m);
	tmp_mx=mx;
	null_mx=k_mx_null_matrix(n,m);
	if(row_num01<0 || row_num01>=n)	return null_mx;
	if(row_num02<0 || row_num02>=n)	return null_mx;

	tmp_mx[row_num01]=mx[row_num02];
	tmp_mx[row_num02]=mx[row_num01];

	return tmp_mx;
}
/************************************************************************/
// returns matrix of which col_num01-th row is exchanged by col_num02-th row
// CAUTION!	row_num and col_num start from "0" not "1"
/************************************************************************/
matrix<double> k_mx_col_exchange(matrix<double> &mx, int col_num01, int col_num02)
{
	int m,n,i;
	n=mx.numrows();
	m=mx.numcols();
	matrix<double> tmp_mx(n,m);
	matrix<double> null_mx(n,m);
	tmp_mx=mx;
	null_mx=k_mx_null_matrix(n,m);
	if(col_num01<0 || col_num01>=m)	return null_mx;
	if(col_num02<0 || col_num02>=m)	return null_mx;

	for(i=0;i<n;i++)
	{
		tmp_mx[i][col_num01]=mx[i][col_num02];
		tmp_mx[i][col_num02]=mx[i][col_num01];
	}
	return tmp_mx;
}
/************************************************************************/
// extract col_num-th column of given matrix
//
/************************************************************************/
vector<double> k_mx_column_extract(matrix<double> &mx, int col_num)
{
	int m,n,i;
	n=mx.numrows();
	m=mx.numcols();
	vector<double> tmp_vt(n);
	for(i=0;i<n;i++)
		tmp_vt[i]=mx[i][col_num];

	return tmp_vt;

}