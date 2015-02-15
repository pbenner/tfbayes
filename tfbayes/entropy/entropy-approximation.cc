/* Copyright (C) 2015 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <tfbayes/utility/probability.hh>

#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/dirichlet_distribution.hpp>
#include <boost/math/distributions/dirichlet.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <tfbayes/utility/histogram.hh>
#include <tfbayes/utility/newton.hh>
#include <tfbayes/utility/summation.hh>
#include <tfbayes/entropy/entropy.hh>

// type declarations
////////////////////////////////////////////////////////////////////////////////
typedef long double real_t;
typedef probability_t<real_t> p_t;
typedef std::vector<p_t> p_vector_t;
typedef histogram_t<real_t, p_t> hist_t;
////////////////////////////////////////////////////////////////////////////////

template <typename T>
void seed_rng(T& rng)
{
        struct timeval tv;
        gettimeofday(&tv, NULL);
        rng.seed(tv.tv_sec*tv.tv_usec);
}

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
        for (size_t i = 0; i < v.size(); i++) {
                if (i != 0) o << " ";
                o << std::fixed << std::setprecision(8) << v[i];
        }
        return o;
}

using namespace std;

void
save_table(const hist_t& histogram, size_t k)
{
        ofstream ofs((boost::format("entropy-approximation-%d.csv") % k).str());

        BOOST_FOREACH(const real_t& x, histogram.x()) {
                ofs << boost::format("%0.8f %e") % (x*std::log(k)) % std::log(pdf(histogram, x))
                    << endl;
        }
}

void
save_histogram(const hist_t& histogram, size_t k)
{
        string filename = (boost::format("entropy-approximation-%d.hh") % k).str();
        ofstream ofs(filename);

        ofs << "/* This file was automatically generated. */"
            << endl << endl;
        ofs << boost::format("std::vector<double> entropy_histogram_%d = {") % k
            << endl;
        for (size_t i = 0; i < histogram.size(); i++) {
                ofs << setprecision(100)
                    << "\t" << std::log(histogram[i]);
                if (i+1 == histogram.size()) {
                        ofs << endl;
                }
                else {
                        ofs << "," << endl;
                }
        }
        ofs << "};" << endl;
        ofs << boost::format("std::vector<double> entropy_histogram_%d_counts = {") % k
            << endl;
        for (size_t i = 0; i < histogram.size(); i++) {
                ofs << setprecision(100)
                    << "\t" << histogram.counts()[i];
                if (i+1 == histogram.size()) {
                        ofs << endl;
                }
                else {
                        ofs << "," << endl;
                }
        }
        ofs << "};" << endl;
}

class
proposal_distribution_t
{
        typedef boost::random::dirichlet_distribution<real_t, p_t> rdirichlet_t;
        typedef boost::math  ::dirichlet_distribution<real_t     > ddirichlet_t;

        // cardinality
        size_t m_k;
        // number of components
        size_t m_size;
        // use exact summation
        bool m_extended_precision;
        // sampling distributions and densities
        std::vector<rdirichlet_t> m_rdirichlet;
        std::vector<ddirichlet_t> m_ddirichlet;

        real_t m_mean(real_t alpha) {
                alpha = std::max(alpha, real_t(1.0e-10));
                return boost::math::digamma(m_k*alpha + 1.0) - boost::math::digamma(alpha + 1.0);
        }
        real_t m_dmean(real_t alpha) {
                alpha = std::max(alpha, real_t(1.0e-10));
                return m_k*boost::math::trigamma(m_k*alpha + 1.0) - boost::math::trigamma(alpha + 1.0);
        }
        real_t m_sigma(real_t alpha) {
                alpha = std::max(alpha, real_t(1.0e-10));
                return std::sqrt((alpha+1.0)/(m_k*alpha+1.0)*boost::math::trigamma(alpha + 1.0)
                                 - boost::math::trigamma(m_k*alpha + 1.0));
        }
        real_t compute_alpha(real_t alpha, real_t target) {
                return newton<real_t>(boost::bind(&proposal_distribution_t::m_mean,  this, _1),
                                      boost::bind(&proposal_distribution_t::m_dmean, this, _1),
                                      alpha, target);
        }
public:
        proposal_distribution_t(size_t k, const hist_t& histogram, bool extended_precision = false, real_t n = 1.0)
                : m_k(k), m_size(0.0), m_extended_precision(extended_precision) {

                real_t alpha_min = std::max(real_t(0.05), compute_alpha(1.0, histogram.x()[0]));
                real_t alpha_max = compute_alpha(1.0, histogram.x()[histogram.size()-1]-histogram.width()/0.6);

                for (real_t alpha = alpha_max; alpha > alpha_min;) {
                        // verbose
                        cout << boost::format("Adding distribution at alpha = %f (with mean entropy %f)")
                                % alpha % m_mean(alpha) << endl;
                        // add distributions
                        m_rdirichlet.push_back(rdirichlet_t(k, alpha));
                        m_ddirichlet.push_back(ddirichlet_t(k, alpha));
                        // compute new alpha
                        alpha = compute_alpha(alpha, m_mean(alpha) - n*m_sigma(alpha));
                        // increase number of components
                        m_size += 1;
                }
                cout << "Done." << endl;
        }

        template <class Engine>
        p_vector_t operator()(Engine& eng) {
                boost::random::uniform_int_distribution<> rint(0, m_size-1);
                return m_rdirichlet[rint(eng)](eng, m_extended_precision);
        }
        const rdirichlet_t& rdirichlet(size_t i) const {
                return m_rdirichlet[i];
        }
        const ddirichlet_t& ddirichlet(size_t i) const {
                return m_ddirichlet[i];
        }
        const size_t& size() const {
                return m_size;
        }
};

p_t pdf(const proposal_distribution_t& d, const p_vector_t& x) {
        p_t result = 0.0;
                
        for (size_t i = 0; i < d.size(); i++) {
                result += from_log_scale(boost::math::log_pdf(d.ddirichlet(i), x));
        }
        return result/p_t(d.size());
}

hist_t
approximate_distribution(size_t k, size_t minimum_counts, size_t bins, bool extended_precision = false)
{
        boost::random::mt19937 gen; seed_rng(gen);
        hist_t histogram(0.0, log(k), bins);
        proposal_distribution_t proposal_distribution(k, histogram, extended_precision);
        p_vector_t theta;
        real_t x;

        for (size_t i = 0; histogram.min_counts() < minimum_counts; i++) {
                while (true) {
                        try {
                                theta = proposal_distribution(gen);
                                break;
                        }
                        catch (std::domain_error &e) {
                                cerr << e.what()
                                     << endl;;
                        }
                }
                x = static_cast<real_t>(entropy(theta, extended_precision))/std::log(k);
                histogram.add(x, 1.0/pdf(proposal_distribution, theta));
                if ((i+1) % 100000 == 0) {
                        vector<real_t>::const_iterator it =
                                std::min_element(histogram.counts().begin(),
                                                 histogram.counts().end());
                        cout << boost::format("-> min counts: %f at %d (%f)")
                                % *it % (it - histogram.counts().begin()) % histogram.x()[it - histogram.counts().begin()]
                             << endl;
                        // save partial results
                        save_table    (histogram, k);
                        save_histogram(histogram, k);
                }
        }
        return histogram;
}

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s from [to]\n\n", pname);
}

int
main(int argc, char *argv[])
{
        if (argc != 2 && argc != 3) {
                print_usage(argv[0], stderr);
                exit(EXIT_FAILURE);
        }
        const size_t from = atoi(argv[1]);
        const size_t to   = argc == 2 ? from : atoi(argv[2]);

        const size_t minimum_samples = 500000;
        const size_t bins = 100;

        for (size_t k = from; k <= to; k++) {
                cerr << boost::format("Sampling entropies for cardinality %d...") % k
                     << endl;

                const hist_t histogram = approximate_distribution(k, minimum_samples, bins);

                save_table    (histogram, k);
                save_histogram(histogram, k);
        }
}
