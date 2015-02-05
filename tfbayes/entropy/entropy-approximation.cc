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

#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <tfbayes/utility/histogram.hh>
#include <tfbayes/utility/probability.hh>
#include <tfbayes/utility/distribution.hh>
#include <tfbayes/entropy/entropy.hh>

typedef histogram_t<double> prob_histogram_t;

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
save_table(const prob_histogram_t& histogram, size_t k)
{
        ofstream ofs((boost::format("entropy-approximation-%d.csv") % k).str());

        BOOST_FOREACH(const double& x, histogram.x()) {
                ofs << boost::format("%0.8f %0.8f") % x % histogram.pdf(x)
                    << endl;
        }
}

void
save_histogram(const prob_histogram_t& histogram, size_t k)
{
        string filename = (boost::format("entropy-approximation-%d.hh") % k).str();
        ofstream ofs(filename);

        ofs << "/* This file was automatically generated. */"
            << endl << endl;
        ofs << boost::format("std::vector<double> entropy_histogram_%d_counts = {") % k
            << endl;
        for (size_t i = 0; i < histogram.size(); i++) {
                ofs << setprecision(100)
                    << "\t" << histogram[i];
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
        size_t m_k;
        double m_lambda;

        boost::random::dirichlet_distribution<double, probability_t> m_rdir1;
        boost::random::dirichlet_distribution<double, probability_t> m_rdir2;
        boost::random::dirichlet_distribution<double, probability_t> m_rdir3;
        boost::math::dirichlet_distribution<> m_ddir1;
        boost::math::dirichlet_distribution<> m_ddir2;
        boost::math::dirichlet_distribution<> m_ddir3;
public:
        proposal_distribution_t(size_t k, double alpha, double beta, double gamma)
                : m_k      (k),
                  m_rdir1   (std::vector<double>(k, alpha)),
                  m_rdir2   (std::vector<double>(k, beta )),
                  m_rdir3   (std::vector<double>(k, gamma)),
                  m_ddir1   (std::vector<double>(k, alpha)),
                  m_ddir2   (std::vector<double>(k, beta )),
                  m_ddir3   (std::vector<double>(k, gamma))
                { }

        template <class Engine>
        std::vector<probability_t> operator()(Engine& eng) {
                boost::random::uniform_int_distribution<> rint(1,3);
                switch(rint(eng)) {
                default:
                case 1: return m_rdir1(eng);
                case 2: return m_rdir2(eng);
                case 3: return m_rdir3(eng);
                }
        }

        double pdf(std::vector<probability_t> x) {
                return 1.0/3.0*boost::math::pdf(m_ddir1, x) +
                       1.0/3.0*boost::math::pdf(m_ddir2, x) +
                       1.0/3.0*boost::math::pdf(m_ddir3, x);
        }
};

struct parameters_t
{
        size_t k;
        double alpha;
        double beta;
        double gamma;
};

prob_histogram_t
approximate_distribution(const parameters_t& parameters, size_t minimum_counts, size_t bins)
{
        boost::random::mt19937 gen; seed_rng(gen);
        proposal_distribution_t proposal_distribution(
                parameters.k, parameters.alpha, parameters.beta, parameters.gamma);
        prob_histogram_t histogram(0.0, log(parameters.k), bins);

        for (size_t i = 0; histogram.min_counts() < minimum_counts; i++) {
                vector<probability_t> theta = proposal_distribution(gen);
                histogram.add(entropy(theta), 1.0/proposal_distribution.pdf(theta));
                if ((i+1) % 100000 == 0) {
                        std::vector<double>::const_iterator it =
                                std::min_element(histogram.counts().begin(),
                                                 histogram.counts().end());
                        cout << boost::format("-> min counts: %f at %d") % *it % (it - histogram.counts().begin())
                             << endl;
                }
        }
        return histogram;
}

int
main(void)
{
        const size_t minimum_samples = 100000;
        const size_t bins = 100;

        parameters_t parameters[] = {
//                { 5, 1.0, 0.2, 5.0 },
//                { 6, 1.0, 0.15, 3.0 },
//                { 7, 0.8, 0.08, 10.0 },
                { 50, 1.0, 0.01, 0.9 },
                { 0, 0, 0, 0 }
        };

        for (size_t i = 0; parameters[i].k; i++) {
                cerr << boost::format("Sampling entropies on simplices of dimension %d...") % parameters[i].k
                     << endl;

                const prob_histogram_t histogram = approximate_distribution(parameters[i], minimum_samples, bins);

                save_table    (histogram, parameters[i].k);
                save_histogram(histogram, parameters[i].k);
        }

}
