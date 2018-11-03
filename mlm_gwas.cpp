#include <cmath>
#include <limits>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>
#include "cmdline.h"
#include "vcf.h"
#include "pheno.h"
#include "statsutil.h"
#include "emma.h"
#include "lsfit.h"
#include "statsutil.h"
#include "lapack.h"
#include "util.h"


using std::ptrdiff_t;
using std::size_t;


namespace {


struct Parameter
{
    std::string vcf;
    std::string pheno;
    std::string covar;
    std::string kin;
    std::string out;
    bool exact = false;
    bool openmp = false;
} par ;


void filter_ind(const std::vector<size_t> &idx, Genotype &gt)
{
    if (gt.ploidy == 1) {
        for (auto &v : gt.dat)
            subset(v,idx).swap(v);
    }
    else {
        std::vector<size_t> idx2;
        for (auto i : idx) {
            idx2.push_back(i*2);
            idx2.push_back(i*2+1);
        }
        for (auto &v : gt.dat)
            subset(v,idx2).swap(v);
    }

    subset(gt.ind, idx).swap(gt.ind);
}

void filter_ind(const std::vector<size_t> &idx, SquareData &sd)
{
    subset(sd.ind,idx).swap(sd.ind);
    subset(sd.dat,idx).swap(sd.dat);
    for (auto &e : sd.dat)
        subset(e,idx).swap(e);
}

void filter_ind(const std::vector<size_t> &idx, Phenotype &pt)
{
    subset(pt.ind,idx).swap(pt.ind);

    if ( ! pt.env.empty() )
        subset(pt.env,idx).swap(pt.env);

    if ( ! pt.blk.empty() )
        subset(pt.blk,idx).swap(pt.blk);

    for (auto &e : pt.dat)
        subset(e,idx).swap(e);
}

void filter_ind(const std::vector<size_t> &idx, Covariate &ct)
{
    subset(ct.ind,idx).swap(ct.ind);

    for (auto &e : ct.dat)
        subset(e,idx).swap(e);
}

void merge(Genotype &gt, SquareData &kin, Phenotype &pt, Covariate &ct, std::vector<size_t> &gi)
{
    bool docovar = ! ct.phe.empty() && ! ct.ind.empty();

    bool domerge = gt.ind != kin.ind || gt.ind != pt.ind;
    if ( ! domerge && docovar )
        domerge = gt.ind != ct.ind;

    if ( ! domerge ) {
        gi.resize(gt.ind.size());
        std::iota(gi.begin(), gi.end(), 0);
        return;
    }

    std::cerr << "INFO: performing data intersection by individual...\n";

    if (gt.ind != kin.ind) {
        std::vector<size_t> idx1, idx2;
        for (auto &e : intersect(gt.ind, kin.ind)) {
            idx1.push_back( index(gt.ind, e) );
            idx2.push_back( index(kin.ind, e) );
        }
        filter_ind(idx1, gt);
        filter_ind(idx2, kin);
    }

    auto ind = intersect(gt.ind, pt.ind);
    if ( docovar )
        ind = intersect(ind, ct.ind);

    gi.clear();
    std::vector<size_t> pi, ci;

    auto n = pt.ind.size();
    for (size_t i = 0; i < n; ++i) {
        if ( std::binary_search(ind.begin(), ind.end(), pt.ind[i]) ) {
            pi.push_back(i);
            gi.push_back( index(gt.ind, pt.ind[i]) );
            if ( docovar )
                ci.push_back( index(ct.ind, pt.ind[i]) );
        }
    }

    filter_ind(pi, pt);

    if ( docovar )
        filter_ind(ci, ct);

    std::cerr << "INFO: there are " << ind.size() << " individuals after intersection\n";
}

void parse_envblk(const Phenotype &pt, std::vector< std::vector<double> > &ac, std::vector< std::vector<double> > &ic)
{
    std::vector< std::vector<double> > xenv, xblk;

    if ( ! pt.env.empty() ) {
        idummy1(factor(pt.env), xenv);
        ac.insert(ac.end(), xenv.begin(), xenv.end());
        ic.insert(ic.end(), xenv.begin(), xenv.end());
    }

    if ( ! pt.blk.empty() ) {
        idummy1(factor(pt.blk), xblk);
        if ( xenv.empty() )
            ac.insert(ac.end(), xblk.begin(), xblk.end());
        else {
            idummy3(factor(pt.env), xenv);
            for (auto &e : xenv) {
                for (auto v : xblk) {
                    std::transform(e.begin(), e.end(), v.begin(), v.begin(), std::multiplies<double>());
                    ac.push_back(v);
                }
            }
        }
    }
}

int assoc_mlm(const Genotype &gt, const SquareData &kin, const std::vector<size_t> &gi,
              const std::vector<double> &y,
              const std::vector< std::vector<double> > &ac,
              std::vector<double> &ps)
{
    auto n = y.size();
    auto m = gt.dat.size();

    ps.assign(m, std::numeric_limits<double>::quiet_NaN());

    auto q0 = 1 + ac.size();
    std::vector<double> x0(n, 1);
    for (auto &v : ac)
        x0.insert(x0.end(), v.begin(), v.end());

    std::vector<double> ki;

    auto K = subset(kin.dat, gi);
    for (auto &e : K) {
        e = subset(e, gi);
        ki.insert(ki.end(), e.begin(), e.end());
    }

    for (size_t i = 0; i < n; ++i)
        ki[i*n+i] += 1;

    EMMA emma;

    int info = emma.solve(n, q0, x0.data(), y.data(), ki.data());

    if (info != 0) {
        std::cerr << "ERROR: failed to solve mixed model\n";
        return 1;
    }

    std::vector<size_t> idx;
    std::vector<allele_t> g1;
    std::vector< std::pair<allele_t,allele_t> > g2;
    std::vector<double> x;
    std::vector<double> v;

    for (size_t j = 0; j < m; ++j) {
        idx.clear();
        std::vector< std::vector<double> > x1;

        if (gt.ploidy == 1) {
            g1.clear();
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii];
                if ( a ) {
                    g1.push_back(a);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g1), x1);
        }
        else {
            g2.clear();
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii*2];
                auto b = gt.dat[j][ii*2+1];
                if ( a && b ) {
                    if (a > b)
                        std::swap(a, b);
                    g2.emplace_back(a, b);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g2), x1);
        }

        if ( x1.empty() )
            continue;

        auto nv = idx.size();
        auto q1 = x1.size();
        auto q = q0 + q1;

        if (nv <= q) {
            std::cerr << "ERROR: not enough observations for " << j + 1 << "\n";
            continue;
        }

        double fval = 0;

        if (nv == n) {
            x = x0;

            for (auto &e : x1)
                x.insert(x.end(), e.begin(), e.end());

            v = ki;
            for (size_t i = 0; i < n; ++i)
                v[i*n+i] -= 1;
            C_dscal(v.size(), emma.vg, v.data(), 1);
            for (size_t i = 0; i < n; ++i)
                v[i*n+i] += emma.ve;

            fval = glsfstat(q1, y, x, v);
        }
        else {
            x.assign(nv, 1);

            for (auto &e : ac) {
                auto z = subset(e, idx);
                x.insert(x.end(), z.begin(), z.end());
            }

            for (auto &e : x1)
                x.insert(x.end(), e.begin(), e.end());

            v.clear();
            for (auto &e : subset(K, idx)) {
                auto z = subset(e, idx);
                v.insert(v.end(), z.begin(), z.end());
            }
            C_dscal(v.size(), emma.vg, v.data(), 1);
            for (size_t i = 0; i < nv; ++i)
                v[i*nv+i] += emma.ve;

            fval = glsfstat(q1, subset(y,idx), x, v);
        }

        ps[j] = fpval(fval, q1, n-q1);
    }

    return 0;
}

int assoc_mlm_exact(const Genotype &gt, const SquareData &kin, const std::vector<size_t> &gi,
                    const std::vector<double> &y,
                    const std::vector< std::vector<double> > &ac,
                    std::vector<double> &ps)
{
    auto n = y.size();
    auto m = gt.dat.size();

    ps.assign(m, std::numeric_limits<double>::quiet_NaN());

    auto q0 = 1 + ac.size();
    std::vector< std::vector<double> > x0;
    x0.push_back( std::vector<double>(n, 1) );
    x0.insert(x0.end(), ac.begin(), ac.end());

    std::vector<double> ki;

    auto K = subset(kin.dat, gi);
    for (auto &e : K) {
        e = subset(e, gi);
        ki.insert(ki.end(), e.begin(), e.end());
    }

    for (size_t i = 0; i < n; ++i)
        ki[i*n+i] += 1;

    std::vector<size_t> idx;
    std::vector<allele_t> g1;
    std::vector< std::pair<allele_t,allele_t> > g2;
    std::vector<double> x;
    std::vector<double> v;
    std::vector<double> ki2;

    for (size_t j = 0; j < m; ++j) {
        idx.clear();
        std::vector< std::vector<double> > x1;

        if (gt.ploidy == 1) {
            g1.clear();
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii];
                if ( a ) {
                    g1.push_back(a);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g1), x1);
        }
        else {
            g2.clear();
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii*2];
                auto b = gt.dat[j][ii*2+1];
                if ( a && b ) {
                    if (a > b)
                        std::swap(a, b);
                    g2.emplace_back(a, b);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g2), x1);
        }

        if ( x1.empty() )
            continue;

        auto nv = idx.size();
        auto q1 = x1.size();
        auto q = q0 + q1;

        if (nv <= q) {
            std::cerr << "ERROR: not enough observations for " << j + 1 << "\n";
            continue;
        }

        x.clear();

        if (nv != n) {
            for (auto &e : x0) {
                auto z = subset(e, idx);
                x.insert(x.end(), z.begin(), z.end());
            }
        }
        else {
            for (auto &e : x0)
                x.insert(x.end(), e.begin(), e.end());
        }

        for (auto &e : x1)
            x.insert(x.end(), e.begin(), e.end());

        EMMA emma;
        double fval = 0;

        if (nv == n) {
            int info = emma.solve(n, q, x.data(), y.data(), ki.data());

            if (info != 0) {
                std::cerr << "ERROR: failed to solve mixed model for " << j + 1 << "\n";
                continue;
            }

            v = ki;
            for (size_t i = 0; i < n; ++i)
                v[i*n+i] -= 1;
            C_dscal(v.size(), emma.vg, v.data(), 1);
            for (size_t i = 0; i < n; ++i)
                v[i*n+i] += emma.ve;

            fval = glsfstat(q1, y, x, v);
        }
        else {
            auto y2 = subset(y, idx);

            ki2.clear();
            for (auto &e : subset(K, idx)) {
                auto z = subset(e, idx);
                ki2.insert(ki2.end(), z.begin(), z.end());
            }
            for (size_t i = 0; i < nv; ++i)
                ki2[i*nv+i] += 1;

            int info = emma.solve(nv, q, x.data(), y2.data(), ki2.data());

            if (info != 0) {
                std::cerr << "ERROR: failed to solve mixed model for " << j + 1 << "\n";
                continue;
            }

            v = ki2;
            for (size_t i = 0; i < nv; ++i)
                v[i*nv+i] -= 1;
            C_dscal(v.size(), emma.vg, v.data(), 1);
            for (size_t i = 0; i < nv; ++i)
                v[i*nv+i] += emma.ve;

            fval = glsfstat(q1, y2, x, v);
        }

        ps[j] = fpval(fval, q1, n-q1);
    }

    return 0;
}

int assoc_mlm_omp(const Genotype &gt, const SquareData &kin, const std::vector<size_t> &gi,
                  const std::vector<double> &y,
                  const std::vector< std::vector<double> > &ac,
                  std::vector<double> &ps)
{
    auto n = y.size();
    auto m = gt.dat.size();

    ps.assign(m, std::numeric_limits<double>::quiet_NaN());

    auto q0 = 1 + ac.size();
    std::vector<double> x0(n, 1);
    for (auto &v : ac)
        x0.insert(x0.end(), v.begin(), v.end());

    std::vector<double> ki;

    auto K = subset(kin.dat, gi);
    for (auto &e : K) {
        e = subset(e, gi);
        ki.insert(ki.end(), e.begin(), e.end());
    }

    for (size_t i = 0; i < n; ++i)
        ki[i*n+i] += 1;

    EMMA emma;

    int info = emma.solve(n, q0, x0.data(), y.data(), ki.data());

    if (info != 0) {
        std::cerr << "ERROR: failed to solve mixed model\n";
        return 1;
    }

    // in earlier OpenMP specifications (<3.0), unsigned integer is not allowed in loop construct
    auto m2 = static_cast<ptrdiff_t>(m);

    #pragma omp parallel for
    for (ptrdiff_t j2 = 0; j2 < m2; ++j2) {
        auto j = static_cast<size_t>(j2);

        std::vector<size_t> idx;
        std::vector< std::vector<double> > x1;

        if (gt.ploidy == 1) {
            std::vector<allele_t> g;
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii];
                if ( a ) {
                    g.push_back(a);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g), x1);
        }
        else {
            std::vector< std::pair<allele_t,allele_t> > g;
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii*2];
                auto b = gt.dat[j][ii*2+1];
                if ( a && b ) {
                    if (a > b)
                        std::swap(a, b);
                    g.emplace_back(a, b);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g), x1);
        }

        if ( x1.empty() )
            continue;

        auto nv = idx.size();
        auto q1 = x1.size();
        auto q = q0 + q1;

        if (nv <= q) {
            std::cerr << "ERROR: not enough observations for " << j + 1 << "\n";
            continue;
        }

        double fval = 0;

        if (nv == n) {
            auto x = x0;

            for (auto &e : x1)
                x.insert(x.end(), e.begin(), e.end());

            auto v = ki;
            for (size_t i = 0; i < n; ++i)
                v[i*n+i] -= 1;
            C_dscal(v.size(), emma.vg, v.data(), 1);
            for (size_t i = 0; i < n; ++i)
                v[i*n+i] += emma.ve;

            fval = glsfstat(q1, y, x, v);
        }
        else {
            std::vector<double> x(nv, 1);

            for (auto &e : ac) {
                auto z = subset(e, idx);
                x.insert(x.end(), z.begin(), z.end());
            }

            for (auto &e : x1)
                x.insert(x.end(), e.begin(), e.end());

            std::vector<double> v;
            for (auto &e : subset(K, idx)) {
                auto z = subset(e, idx);
                v.insert(v.end(), z.begin(), z.end());
            }
            C_dscal(v.size(), emma.vg, v.data(), 1);
            for (size_t i = 0; i < nv; ++i)
                v[i*nv+i] += emma.ve;

            fval = glsfstat(q1, subset(y, idx), x, v);
        }

        ps[j] = fpval(fval, q1, n-q1);
    }

    return 0;
}

int assoc_mlm_exact_omp(const Genotype &gt, const SquareData &kin, const std::vector<size_t> &gi,
                        const std::vector<double> &y,
                        const std::vector< std::vector<double> > &ac,
                        std::vector<double> &ps)
{
    auto n = y.size();
    auto m = gt.dat.size();

    ps.assign(m, std::numeric_limits<double>::quiet_NaN());

    auto q0 = 1 + ac.size();
    std::vector< std::vector<double> > x0;
    x0.push_back( std::vector<double>(n, 1) );
    x0.insert(x0.end(), ac.begin(), ac.end());

    std::vector<double> ki;

    auto K = subset(kin.dat, gi);
    for (auto &e : K) {
        e = subset(e, gi);
        ki.insert(ki.end(), e.begin(), e.end());
    }

    for (size_t i = 0; i < n; ++i)
        ki[i*n+i] += 1;

    // in earlier OpenMP specifications (<3.0), unsigned integer is not allowed in loop construct
    auto m2 = static_cast<ptrdiff_t>(m);

    #pragma omp parallel for
    for (ptrdiff_t j2 = 0; j2 < m2; ++j2) {
        auto j = static_cast<size_t>(j2);

        std::vector<size_t> idx;
        std::vector< std::vector<double> > x1;

        if (gt.ploidy == 1) {
            std::vector<allele_t> g;
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii];
                if ( a ) {
                    g.push_back(a);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g), x1);
        }
        else {
            std::vector< std::pair<allele_t,allele_t> > g;
            for (size_t i = 0; i < n; ++i) {
                auto ii = gi[i];
                auto a = gt.dat[j][ii*2];
                auto b = gt.dat[j][ii*2+1];
                if ( a && b ) {
                    if (a > b)
                        std::swap(a, b);
                    g.emplace_back(a, b);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g), x1);
        }

        if ( x1.empty() )
            continue;

        auto nv = idx.size();
        auto q1 = x1.size();
        auto q = q0 + q1;

        if (nv <= q) {
            std::cerr << "ERROR: not enough observations for " << j + 1 << "\n";
            continue;
        }

        std::vector<double> x;

        if (nv != n) {
            for (auto &e : x0) {
                auto z = subset(e, idx);
                x.insert(x.end(), z.begin(), z.end());
            }
        }
        else {
            for (auto &e : x0)
                x.insert(x.end(), e.begin(), e.end());
        }

        for (auto &e : x1)
            x.insert(x.end(), e.begin(), e.end());

        EMMA emma;
        double fval = 0;

        if (nv == n) {
            int info = emma.solve(n, q, x.data(), y.data(), ki.data());

            if (info != 0) {
                std::cerr << "ERROR: failed to solve mixed model for " << j + 1 << "\n";
                continue;
            }

            auto v = ki;
            for (size_t i = 0; i < n; ++i)
                v[i*n+i] -= 1;
            C_dscal(v.size(), emma.vg, v.data(), 1);
            for (size_t i = 0; i < n; ++i)
                v[i*n+i] += emma.ve;

            fval = glsfstat(q1, y, x, v);
        }
        else {
            auto y2 = subset(y, idx);

            std::vector<double> ki2;
            for (auto &e : subset(K, idx)) {
                auto z = subset(e, idx);
                ki2.insert(ki2.end(), z.begin(), z.end());
            }
            for (size_t i = 0; i < nv; ++i)
                ki2[i*nv+i] += 1;

            int info = emma.solve(nv, q, x.data(), y2.data(), ki2.data());

            if (info != 0) {
                std::cerr << "ERROR: failed to solve mixed model for " << j + 1 << "\n";
                continue;
            }

            auto v = ki2;
            for (size_t i = 0; i < nv; ++i)
                v[i*nv+i] -= 1;
            C_dscal(v.size(), emma.vg, v.data(), 1);
            for (size_t i = 0; i < nv; ++i)
                v[i*nv+i] += emma.ve;

            fval = glsfstat(q1, y2, x, v);
        }

        ps[j] = fpval(fval, q1, n-q1);
    }

    return 0;
}


} // namespace


int mlm_gwas(int argc, char *argv[])
{
    std::cerr << "MLM-GWAS 1.0 (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd;

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--pheno", "phenotype file", "");
    cmd.add("--covar", "covariate file", "");
    cmd.add("--kin", "kinship file", "");
    cmd.add("--out", "output file", "mlm-gwas.out");
    cmd.add("--exact", "using exact test statistic");
    cmd.add("--openmp", "enable OpenMP multithreading");

    cmd.parse(argc, argv);

    if (argc < 2) {
        cmd.show();
        return 1;
    }

    par.vcf = cmd.get("--vcf");
    par.pheno = cmd.get("--pheno");
    par.covar = cmd.get("--covar");
    par.kin = cmd.get("--kin");
    par.out = cmd.get("--out");
    par.exact = cmd.has("--exact");
    par.openmp = cmd.has("--openmp");

    Genotype gt;
    Phenotype pt;
    Covariate ct;
    SquareData kin;

    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    std::cerr << "INFO: reading phenotype file...\n";
    if (read_pheno(par.pheno, pt) != 0)
        return 1;
    std::cerr << "INFO: " << pt.ind.size() << " observations, " << pt.phe.size() << " traits\n";

    std::cerr << "INFO: reading kinship file...\n";
    if (read_square(par.kin, kin) != 0)
        return 1;
    std::cerr << "INFO: " << kin.ind.size() << " individuals\n";

    if ( ! par.covar.empty() ) {
        std::cerr << "INFO: reading covariate file...\n";
        if (read_covar(par.covar, ct) != 0)
            return 1;
        std::cerr << "INFO: " << ct.ind.size() << " individuals, " << ct.phe.size() << " covariates\n";
    }

    std::vector<size_t> gi;
    merge(gt, kin, pt, ct, gi);

    if ( gi.empty() ) {
        std::cerr << "ERROR: no valid observations are found\n";
        return 1;
    }

    std::vector< std::vector<double> > ac, ic;
    parse_envblk(pt, ac, ic);

    ac.insert(ac.end(), ct.dat.begin(), ct.dat.end());

    std::vector< std::vector<double> > vps;

    auto nt = pt.phe.size();

    for (size_t t = 0; t < nt; ++t) {
        std::vector<char> keep;
        for (auto e : pt.dat[t])
            keep.push_back( std::isfinite(e) );

        auto y = pt.dat[t];
        auto ac2 = ac;
        auto gi2 = gi;

        if (std::find(keep.begin(), keep.end(), 0) != keep.end()) {
            y = subset(y, keep);
            for (auto &e : ac2)
                e = subset(e, keep);
            gi2 = subset(gi, keep);
        }

        std::vector<double> ps;

        if ( par.exact ) {
            if ( par.openmp )
                assoc_mlm_exact_omp(gt, kin, gi2, y, ac2, ps);
            else
                assoc_mlm_exact(gt, kin, gi2, y, ac2, ps);
        }
        else {
            if ( par.openmp )
                assoc_mlm_omp(gt, kin, gi2, y, ac2, ps);
            else
                assoc_mlm(gt, kin, gi2, y, ac2, ps);
        }

        vps.push_back(ps);
    }

    std::ofstream ofs(par.out + ".pval");

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << par.out << ".pval" << "\n";
        return 1;
    }

    ofs << "Locus\tChromosome\tPosition";
    for (auto &e : pt.phe)
        ofs << "\t" << e;
    ofs << "\n";

    auto m = gt.loc.size();

    for (size_t j = 0; j < m; ++j) {
        ofs << gt.loc[j] << "\t" << gt.chr[j] << "\t" << gt.pos[j];
        for (size_t t = 0; t < nt; ++t)
            ofs << "\t" << vps[t][j];
        ofs << "\n";
    }

    return 0;
}
