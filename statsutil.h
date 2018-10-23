#ifndef STATSUTIL_H
#define STATSUTIL_H

#include <vector>
#include <algorithm>


// upper tail F cdf
double fpval(double x, double df1, double df2);


// Create dummy variables for grouping variable g (integers: 0,1,...)
//
// 0/-1/1 coding, full rank, drop last level
int idummy1(const std::vector<int> &g, std::vector< std::vector<double> > &x);

// 0/1 coding, full rank, drop last level
int idummy2(const std::vector<int> &g, std::vector< std::vector<double> > &x);

// 0/1 coding, overdetermined
int idummy3(const std::vector<int> &g, std::vector< std::vector<double> > &x);


// generate factor variable
template<typename T>
std::vector<int> factor(const std::vector<T> &v)
{
    auto u = v;
    std::sort(u.begin(), u.end());
    u.erase(std::unique(u.begin(), u.end()), u.end());

    std::vector<int> gi;

    gi.reserve(v.size());
    for (auto &e : v) {
        auto itr = std::find(u.begin(), u.end(), e);
        gi.push_back(itr - u.begin());
    }

    return gi;
}


// generate factor variable with level names
template<typename T>
void factor(const std::vector<T> &v, std::vector<T> &gn, std::vector<int> &gi)
{
    auto u = v;
    std::sort(u.begin(), u.end());
    u.erase(std::unique(u.begin(), u.end()), u.end());

    gi.clear();
    gi.reserve(v.size());

    for (auto &e : v) {
        auto itr = std::find(u.begin(), u.end(), e);
        gi.push_back(itr - u.begin());
    }

    gn.swap(u);
}


// corrected sum of squares
template<typename T>
double calc_css(const std::vector<T> &x)
{
    double sum = 0;
    for (auto e : x)
        sum += e;
    double mean = sum / x.size();
    double css = 0;
    for (auto e : x)
        css += (e-mean)*(e-mean);
    return css;
}


#endif // STATSUTIL_H
