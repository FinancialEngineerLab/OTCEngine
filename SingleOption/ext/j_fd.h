#ifndef J_FD_H
#define J_FD_H

void gridcontrol(double *px, double *dpx, int minnode, int maxnode, double target, int boundaryflag);
int findlowerindex(double *px, double target, int minnode, int maxnode);
int findnearestindex(const double *px, double target, int minnode, int maxnode);
void trimatrix1d(double *A, double *B, double *C,double *alpha, double *beta,double rfrate,double dt, double *px,double *dpx, int minnode, int maxnode);
void trimxsolve1d(double *A, double *B, double *C, double *vold, double *vnew,int min, int max, int lb, int ub);
int findnearestindex(const double *px, double target, int minnode, int maxnode);

#endif