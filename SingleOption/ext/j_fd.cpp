#include "j_fd.h"
#include <assert.h>
#include <stdexcept>
#include <algorithm>

void trimatrix1d(double *A, double *B, double *C,double *alpha, double *beta,double rfrate,double dt, double *px,double *dpx, int minnode, int maxnode)
{
	for(int i=minnode;i<=maxnode;i++)
	{
       A[i] = -2.0 * alpha[i] * px[i] * px[i] / (dpx[i-1] * (dpx[i] + dpx[i-1])) + beta[i] * px[i] * dpx[i] / (dpx[i-1] * (dpx[i] + dpx[i-1]));
       B[i] = 1.0 + 2.0 * alpha[i] * px[i] * px[i] / (dpx[i] * dpx[i-1]) - beta[i] * px[i] * (dpx[i] - dpx[i-1]) / (dpx[i] * dpx[i-1]) + dt * rfrate;
       C[i] = -2.0 * alpha[i] * px[i] * px[i] / (dpx[i] * (dpx[i] + dpx[i-1])) - beta[i] * px[i] * dpx[i-1] / (dpx[i] * (dpx[i] + dpx[i-1]))	;
	}
}

void trimxsolve1d(double *A, double *B, double *C, double *vold, double *vnew,int min, int max, int lb, int ub)
{
		double* aa=new double[max+1];
		double* bb=new double[max+1];
		double* cc=new double[max+1];
		double* D=new double[max+1];
		double* W=new double[max+1];
		for(int i=min+1;i<=max-1;i++){
            aa[i] = A[i];
            bb[i] = B[i];
            cc[i] = C[i];
            vnew[i] = vold[i];
		}

        //경계조건 처리(for vold)
        //lower boundary
        bb[min + 1] +=2 * (1 - lb) * aa[min + 1];
        cc[min + 1] -= (1 - lb) * aa[min + 1];
        vnew[min + 1] -= lb * aa[min + 1] * vnew[min];

        //upper boundary
        aa[max - 1] -= (1 - ub) * cc[max - 1];
        bb[max - 1] +=2 * (1 - ub) * cc[max - 1];
        vnew[max - 1] -= ub * cc[max - 1] * vnew[max];

        //thomas algorithm
        D[min + 1] = bb[min + 1];
        W[min + 1] = vnew[min + 1];

        for(int i=min+2;i<=max-1;i++){
                D[i] = bb[i] - aa[i] * cc[i - 1] / D[i - 1];
                W[i] = vnew[i] - aa[i] / D[i - 1] * W[i - 1];
		}

        vnew[max - 1] = W[max - 1] / D[max - 1];
        for(int i=max-2;i>=min+1;i--){
			vnew[i] = (W[i] - cc[i] * vnew[i + 1]) / D[i];
		}

        //경계조건 처리(for vnew)
        //upper boundary
        vnew[max] = ub * vnew[max] + (1 - ub) * (2 * vnew[max - 1] - vnew[max - 2]);

        //lower boundary
        vnew[min] = lb * vnew[min] + (1 - lb) * (2 * vnew[min + 1] - vnew[min + 2]);

		delete[] aa;
		delete[] bb;
		delete[] cc;
		delete[] D;
		delete[] W;
}

int findlowerindex(double *px, double target, int minnode, int maxnode)
{//px는 오름차순인 경우

	if(target<=px[minnode])
		return minnode;
	if(target>=px[maxnode])
		return maxnode;

	for(int i=minnode;i<maxnode;i++)
		if(px[i]<=target && target <px[i+1])
			return i;

	assert(0);
	return 0;
}

int findnearestindex(const double *px, double target, int minnode, int maxnode)
{//px는 오름차순정렬
	if(target<=px[minnode])
		return minnode;
	if(target>=px[maxnode])
		return maxnode;
	for(int i=minnode;i<maxnode;i++){
		if(px[i]<=target && target <px[i+1]){
			if(target-px[i] <px[i+1]-target){
				return i;
			}else{
				return i+1;
			}
		}
	}
	throw std::logic_error("interpolaton fail :findnearestindex");
	return -1;

}

void gridcontrol(double *px, double *dpx, int minnode, int maxnode, double target, int boundaryflag)
{
	for(int i=minnode;i<=maxnode;i++)
	{
		if(target >=px[i] && target <=px[i+1])
		{
			if(boundaryflag==0){ //strike
				double tmp=std::min((target-px[i-1])/2, (px[i+2]-target)/2);
				px[i]=target-tmp;
				px[i+1]=target+tmp;
			}else if(boundaryflag==1){  //barrier
				if(target-px[i] <=dpx[i]/2.0)
					px[i]=target;
				else
					px[i+1]=target;
			}

			dpx[i+1]=px[i+2]-px[i+1];
			dpx[i]=px[i+1]-px[i];
			dpx[i-1]=px[i]-px[i-1];
			break;
		}
	}
	//if here, error statement
}