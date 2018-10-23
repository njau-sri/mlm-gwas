#include <cmath>
#include <limits>
#include <memory>
#include <algorithm>
#include "emma.h"
#include "lapack.h"


using std::size_t;


#ifndef M_2PI
#define M_2PI 6.283185307179586476925286766559 // 2*pi
#endif


namespace {

// S   = I - X * (X'X)^-1 * X'
// SHS = S * (K + I) * S
// eigen SHS
int eigen_shs(size_t n, size_t q, const double *x, const double *ki, double *eval, double *evec)
{
    auto qq = q * q;
    auto nq = n * q;
    auto nn = n * n;

    auto lwork = qq + nq + nn*3;
    std::unique_ptr<double[]> work(new double[lwork]);

    auto xx  = work.get();
    auto xxx = xx + qq;
    auto s   = xxx + nq;
    auto sh  = s + nn;
    auto shs = sh + nn;
    auto w   = work.get();
    auto z   = w + n;

    // xx = X'X
    C_dsyrk('U', 'T', q, n, 1.0, x, n, 0.0, xx, q);

    // xx = inv(X'X)
    int info = M_dsyinv(q, xx);
    if (info != 0)
        return 1;

    // xixtx = X*inv(X'X)
    C_dgemm('N', 'N', n, q, q, 1.0, x, n, xx, q, 0.0, xxx, n);

    // S = I
    std::fill_n(s, nn, 0.0);
    for (size_t i = 0; i < n; ++i)
        s[i*n+i] = 1.0;

    // S = I - X*inv(X'X)*X'
    C_dgemm('N', 'T', n, n, q, -1.0, xxx, n, x, n, 1.0, s, n);

    // sk = S*(K+I)
    C_dgemm('N', 'N', n, n, n, 1.0, s, n, ki, n, 0.0, sh, n);

    // sks = S*(K+I)*S
    C_dgemm('N', 'N', n, n, n, 1.0, sh, n, s, n, 0.0, shs, n);

    // eigen S*(K+I)*S
    bint m = 0;
    std::unique_ptr<bint[]> sup(new bint[2*n]);

    info = C_dsyevr('V', 'A', 'U', n, shs, n, 0.0, 0.0, 0, 0, 0.0, &m, w, z, n, sup.get());
    if (info != 0)
        return 2;

    auto p = n - q;
    for (size_t j = 0; j < p; ++j) {
        auto k = n - j - 1;
        eval[j] = w[k] - 1.0;
        for (size_t i = 0; i < n; ++i)
            evec[j*n+i] = z[k*n+i];
    }

    return 0;
}

// REstricted Maximum-Likelihood
//
//   1              n-q                       n-q      eta_s^2         n-q
//   - [ (n-q) LOG ------ - (n-q) - (n-q) LOG SUM ------------------ - SUM LOG(lambda_s + delta) ]
//   2              2*pi                      s=1  lambda_s + delta    s=1
//
// m = n - q
//

double calc_LL_REML(int m, double ldelta, const double *lambda, const double *eta)
{
    double a = 0.0, b = 0.0;
    double delta = std::exp(ldelta);

    for (int i = 0; i < m; ++i) {
        a += eta[i] * eta[i] / (lambda[i] + delta);
        b += std::log(lambda[i] + delta);
    }

    return 0.5 * (m * (std::log(m / M_2PI) - 1 - std::log(a)) - b);
}

// Derivative of REstricted Maximum-Likelihood
//
//    n-q     SUM_s eta_s^2 / (lambda_s + delta)^2    1            1
//   ----- * -------------------------------------- - - SUM ------------------
//     2       SUM_s eta_s^2 / (lambda_s + delta)     2  s   lambda_s + delta
//
// m = n - q
//

double calc_dLL_REML(int m, double ldelta, const double *lambda, const double *eta)
{
    double delta = std::exp(ldelta);
    double a = 0.0, b = 0.0, c = 0.0;

    for (int i = 0; i < m; ++i) {
        double eta2 = eta[i] * eta[i];
        double lamdel = lambda[i] + delta;
        a += eta2 / (lamdel * lamdel);
        b += eta2 / lamdel;
        c += 1.0 / lamdel;
    }

    return 0.5 * (m * (a / b) - c);
}

double bisect_root_REML(double a, double b, double fa, double fb, double tol, int maxit,
                        int m, const double *lambda, const double *eta)
{
    if (fa == 0.0)
        return a;

    if (fb == 0.0)
        return b;

    if (a >= b || fa * fb > 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    for (int i = 0; i < maxit; ++i) {
        if (std::fabs(a - b) < tol)
            break;

        double x = (a + b) / 2;
        double fx = calc_dLL_REML(m, x, lambda, eta);

        if (x == b || x == a)
            break;

        if (fx == 0.0) {
            a = b = x;
            break;
        }
        else if (fx * fa < 0.0) {
            b = x;
            fb = fx;
        }
        else {
            a = x;
            fa = fx;
        }
    }

    return a;
}

} // namespace


int EMMA::solve(int n, int q, const double *x, const double *y, const double *ki)
{
    REML = delta = vg = ve = std::numeric_limits<double>::quiet_NaN();

    if (n <= q)
        return 1;

    int p = n - q;
    int m = grid + 1;

    auto lwork = p + n*p + p + m*3;
    std::unique_ptr<double[]> work(new double[lwork]);

    auto lambda = work.get();
    auto U      = lambda + p;
    auto eta    = U + n*p;
    auto ldelta = eta + p;
    auto LL     = ldelta + m;
    auto dLL    = LL + m;

    int info = eigen_shs(n, q, x, ki, lambda, U);
    if (info != 0)
        return 2;

    // eta = U^T * y
    C_dgemv('T', n, p, 1.0, U, n, y, 1, 0.0, eta, 1);

    for (int i = 0; i < m; ++i) {
        ldelta[i] = llim + i * (ulim - llim) / grid;
        LL[i] = calc_LL_REML(p, ldelta[i], lambda, eta);
        dLL[i] = calc_dLL_REML(p, ldelta[i], lambda, eta);
        if (std::isnan(REML) || LL[i] > REML) {
            REML = LL[i];
            delta = ldelta[i];
        }
    }

    for (int i = 0; i < grid; ++i) {
        if ( dLL[i] > 0.0 && dLL[i+1] < 0.0 && ! std::isnan(LL[i]) ) {
            double root = bisect_root_REML(ldelta[i], ldelta[i+1], dLL[i], dLL[i+1], tol, maxit, p, lambda, eta);
            double optLL = calc_LL_REML(p, root, lambda, eta);
            if ( std::isnan(REML) || optLL > REML ) {
                REML = optLL;
                delta = root;
            }
        }
    }

    if ( ! std::isnan(REML) && ! std::isnan(delta) ) {
        delta = std::exp(delta);
        vg = 0.0;
        for (int i = 0; i < p; ++i)
            vg += eta[i] * eta[i] / (lambda[i] + delta);
        vg /= p;
        ve = vg * delta;
    }

    return 0;
}
