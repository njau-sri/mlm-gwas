#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include "pheno.h"
#include "split.h"


using std::size_t;


int read_pheno(const std::string &filename, Phenotype &pt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    size_t ln = 0;
    std::vector<std::string> colnames;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        split(line, " \t", colnames);

        if ( ! colnames.empty() )
            break;
    }

    std::vector<size_t> jphe;
    size_t jenv = 0, jblk = 0;

    for (size_t j = 1; j < colnames.size(); ++j) {
        if (colnames[j] == "_ENV_")
            jenv = j;
        else if (colnames[j] == "_BLK_")
            jblk = j;
        else
            jphe.push_back(j);
    }

    for (auto j : jphe)
        pt.phe.push_back(colnames[j]);

    std::vector< std::vector<double> > dat;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if ( vs.empty() )
            continue;

        if (vs.size() != colnames.size()) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vs.size()
                      << "!=" << colnames.size() << "\n";
            return 1;
        }

        pt.ind.push_back(vs[0]);

        if (jenv > 0)
            pt.env.push_back(vs[jenv]);

        if (jblk > 0)
            pt.blk.push_back(vs[jblk]);

        std::vector<double> v;

        for (auto j : jphe) {
            if (vs[j] == "?" || vs[j] == "NA" || vs[j] == ".")
                v.push_back(std::numeric_limits<double>::quiet_NaN());
            else
                v.push_back(std::stod(vs[j]));
        }

        dat.push_back(v);
    }

    auto m = pt.phe.size();
    auto n = pt.ind.size();

    for (size_t j = 0; j < m; ++j) {
        std::vector<double> v(n);
        for (size_t i = 0; i < n; ++i)
            v[i] = dat[i][j];
        pt.dat.push_back(v);
    }

    return 0;
}

int write_pheno(const Phenotype &pt, const std::string &filename)
{
    std::ofstream ofs(filename);
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto m = pt.phe.size();
    auto n = pt.ind.size();

    ofs << "Indiv";
    if ( ! pt.env.empty() )
        ofs << "\t_ENV_";
    if ( ! pt.blk.empty() )
        ofs << "\t_BLK_";
    for (size_t j = 0; j < m; ++j)
        ofs << "\t" << pt.phe[j];
    ofs << "\n";

    for (size_t i = 0; i < n; ++i) {
        ofs << pt.ind[i];
        if ( ! pt.env.empty() )
            ofs << "\t" << pt.env[i];
        if ( ! pt.blk.empty() )
            ofs << "\t" << pt.blk[i];
        for (size_t j = 0; j < m; ++j)
            ofs << "\t" << pt.dat[j][i];
        ofs << "\n";
    }

    return 0;
}

int read_covar(const std::string &filename, Covariate &ct)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    size_t ln = 0;
    std::vector<std::string> colnames;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        split(line, " \t", colnames);

        if ( ! colnames.empty() )
            break;
    }

    if (colnames.size() < 2)
        return 0;

    ct.phe.assign(colnames.begin() + 1, colnames.end());

    std::vector< std::vector<double> > dat;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if ( vs.empty() )
            continue;

        if (vs.size() != colnames.size()) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vs.size()
                      << "!=" << colnames.size() << "\n";
            return 1;
        }

        ct.ind.push_back(vs[0]);
        vs.erase(vs.begin());

        std::vector<double> v;

        for (auto &e : vs) {
            v.push_back(std::stod(e));
            if ( ! std::isfinite(v.back()) ) {
                std::cerr << "ERROR: only finite value is allowed: " << e << "\n";
                return 1;
            }
        }

        dat.push_back(v);
    }

    auto m = ct.phe.size();
    auto n = ct.ind.size();

    for (size_t j = 0; j < m; ++j) {
        std::vector<double> v(n);
        for (size_t i = 0; i < n; ++i)
            v[i] = dat[i][j];
        ct.dat.push_back(v);
    }

    return 0;
}

int write_covar(const Covariate &ct, const std::string &filename)
{
    std::ofstream ofs(filename);
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto m = ct.phe.size();
    auto n = ct.ind.size();

    ofs << "Indiv";
    for (size_t j = 0; j < m; ++j)
        ofs << "\t" << ct.phe[j];
    ofs << "\n";

    for (size_t i = 0; i < n; ++i) {
        ofs << ct.ind[i];
        for (size_t j = 0; j < m; ++j)
            ofs << "\t" << ct.dat[j][i];
        ofs << "\n";
    }

    return 0;
}

int read_square(const std::string &filename, SquareData &sd)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    size_t ln = 0, cc = 0;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;

        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if ( vs.empty() )
            continue;

        if (cc == 0)
            cc = vs.size();

        if (vs.size() != cc) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " ("
                      << vs.size() << "!=" << cc << "): " << filename << "\n";
            return 1;
        }

        sd.ind.push_back(vs[0]);
        vs.erase(vs.begin());

        std::vector<double> v;

        for (auto &e : vs) {
            v.push_back(std::stod(e));
            if ( ! std::isfinite(v.back()) ) {
                std::cerr << "ERROR: only finite value is allowed: " << e << "\n";
                return 1;
            }
        }

        sd.dat.push_back(v);
    }

    if (cc != 0 && cc != sd.ind.size() + 1) {
        std::cerr << "ERROR: data must be square: " << filename << "\n";
        return 1;
    }

    return 0;
}

int write_square(const SquareData &sd, const std::string &filename)
{
    std::ofstream ofs(filename);

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto n = sd.ind.size();

    for (size_t i = 0; i < n; ++i) {
        ofs << sd.ind[i];
        for (size_t j = 0; j < n; ++j)
            ofs << "\t" << sd.dat[i][j];
        ofs << "\n";
    }

    return 0;
}
