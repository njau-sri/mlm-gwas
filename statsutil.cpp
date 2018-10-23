#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
#include <boost/math/distributions/fisher_f.hpp>

#include "statsutil.h"


using std::size_t;


double fpval(double x, double df1, double df2)
{
    boost::math::fisher_f f(df1, df2);
    return boost::math::cdf(boost::math::complement(f, x));
}

int idummy1(const std::vector<int> &g, std::vector< std::vector<double> > &x)
{
    if ( g.empty() ) {
        x.clear();
        return 0;
    }

    auto p = std::minmax_element(g.begin(), g.end());

    if (*p.first != 0) {
        x.clear();
        return 1;
    }

    if (*p.second == 0) {
        x.clear();
        return 0;
    }

    auto n = g.size();
    auto m = static_cast<size_t>( * p.second );

    x.assign(m, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        auto k = static_cast<size_t>(g[i]);
        if (k != m)
            x[k][i] = 1.0;
        else
            for (size_t j = 0; j < m; ++j)
                x[j][i] = -1.0;
    }

    return 0;
}

int idummy2(const std::vector<int> &g, std::vector< std::vector<double> > &x)
{
    if ( g.empty() ) {
        x.clear();
        return 1;
    }

    auto p = std::minmax_element(g.begin(), g.end());

    if (*p.first != 0) {
        x.clear();
        return 1;
    }

    if (*p.second == 0) {
        x.clear();
        return 0;
    }

    auto n = g.size();
    auto m = static_cast<size_t>( * p.second );

    x.assign(m, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        auto k = static_cast<size_t>(g[i]);
        if (k != m)
            x[k][i] = 1.0;
    }

    return 0;
}

int idummy3(const std::vector<int> &g, std::vector< std::vector<double> > &x)
{
    if ( g.empty() ) {
        x.clear();
        return 0;
    }

    auto p = std::minmax_element(g.begin(), g.end());

    if (*p.first != 0) {
        x.clear();
        return 1;
    }

    if (*p.second == 0) {
        x.assign(1, std::vector<double>(g.size(), 1.0));
        return 0;
    }

    auto n = g.size();
    auto m = static_cast<size_t>( * p.second );

    x.assign(m + 1, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        auto k = static_cast<size_t>(g[i]);
        x[k][i] = 1.0;
    }

    return 0;
}
