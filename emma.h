#ifndef EMMA_H
#define EMMA_H


// Kang, H.M. et al. Efficient control of population structure in model organism association mapping.
//   Genetics 178, 1709-23 (2008). doi: 10.1534/genetics.107.080101


// EMMA Linear Mixed Model Solver
//
//   y = Xb + Zu + e
//
//   V = vg*ZKZ' + ve*I
//
//   Input
//     n   the number of rows of the matrix X
//     q   the number of columns of the matrix X
//     x   the n*q matrix X
//     y   the n*1 vector y
//     ki  the n*n symmetric matrix K + I
//

struct EMMA
{
    double REML = 0.0;
    double delta = 0.0;
    double vg = 0.0;
    double ve = 0.0;
    double tol = 1e-4;
    double llim = -10.0;
    double ulim = 10.0;
    int grid = 100;
    int maxit = 100;

    int solve(int n, int q, const double *x, const double *y, const double *ki);
};


#endif // EMMA_H
