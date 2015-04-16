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
#include <tfbayes/entropy/entropy-approximation-mixture.hh>
#include <tfbayes/entropy/entropy-approximation-recursive.hh>

// type declarations
////////////////////////////////////////////////////////////////////////////////
typedef long double real_t;
typedef probability_t<real_t> p_t;
typedef std::vector<p_t> p_vector_t;
typedef histogram_t<real_t, p_t> hist_t;
////////////////////////////////////////////////////////////////////////////////

struct parameters_t {
        size_t bins;
        size_t minimum_counts;
        bool extended_precision;
        proposal_distribution_mixture_t<real_t, p_t>::parameters_t mixture_parameters;
};

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

template <typename T, typename Engine>
p_vector_t
draw_proposal(T& proposal_distribution, Engine& gen)
{
        while (true) {
                try {
                        return proposal_distribution(gen);
                }
                catch (std::domain_error &e) {
                        cerr << e.what()
                             << endl;
                }
        }
}

hist_t
approximate_distribution(size_t k, const parameters_t& parameters)
{
        boost::random::mt19937 gen; seed_rng(gen);
        hist_t histogram(0.0, 1.0, parameters.bins);
        // proposal_distribution_mixture_t<real_t, p_t> proposal_distribution(
        //         k,
        //         parameters.mixture_parameters,
        //         parameters.extended_precision);
        proposal_distribution_recursive_t<real_t, p_t> proposal_distribution(k);
        p_vector_t theta;
        real_t x;

        for (size_t i = 0; histogram.min_counts() < parameters.minimum_counts; i++) {
                theta = draw_proposal(proposal_distribution, gen);
                x     = static_cast<real_t>(entropy(theta, parameters.extended_precision))/std::log(k);
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

map<size_t, parameters_t>
default_parameters()
{
        map<size_t, parameters_t> m;

        // alpha_min, alpha_max, sigma_n, bins, samples, precision
        m[ 2] = {100, 500000, false, {1.0, 1.0, 1.0}};
        m[ 3] = {100, 500000, false, {1.0, 1.0, 1.0}};
        m[ 4] = {200, 500000, false, {0.1, 1.5, 1.0}};
        m[ 5] = {200, 500000, false, {0.1, 3.0, 1.0}};
        m[11] = {200, 500000, false, {0.1, 3.0, 1.0}};

        return m;
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

        map<size_t, parameters_t> parameters = default_parameters();

        for (size_t k = from; k <= to; k++) {
                if (parameters.find(k) == parameters.end()) {
                        cerr << boost::format("Error: no parameters found for cardinality %d.")
                                % k
                             << endl;
                        exit(EXIT_FAILURE);
                }
                cerr << boost::format("Sampling entropies for cardinality %d...") % k
                     << endl;

                const hist_t histogram = approximate_distribution(k, parameters[k]);

                save_table    (histogram, k);
                save_histogram(histogram, k);
        }
}
