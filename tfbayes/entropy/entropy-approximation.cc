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
#include <tfbayes/utility/polygamma.hh>
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
        typedef boost::random::dirichlet_distribution<double, probability_t> rdirichlet_t;
        typedef boost::math::dirichlet_distribution<> ddirichlet_t;

        // cardinality
        size_t m_k;
        // number of components
        size_t m_size;
        // sampling distributions and densities
        std::vector<rdirichlet_t> m_rdirichlet;
        std::vector<ddirichlet_t> m_ddirichlet;

        double m_mean(double alpha) {
                return boost::math::digamma(m_k*alpha + 1.0) - boost::math::digamma(alpha + 1.0);
        }
        double m_dmean(double alpha) {
                return m_k*boost::math::trigamma(m_k*alpha + 1.0) - boost::math::trigamma(alpha + 1.0);
        }
        double m_sigma(double alpha) {
                return std::sqrt((alpha+1.0)/(m_k*alpha+1.0)*boost::math::trigamma(alpha + 1.0)
                                 - boost::math::trigamma(m_k*alpha + 1.0));
        }
public:
        proposal_distribution_t(size_t k, double n = 0.5, double alpha_max = 10.0)
                : m_k(k), m_size(0.0) {
                // go to lower alpha values
                for (double alpha = 1.0; alpha > 0.0;) {
                        cout << boost::format("adding distribution at alpha = %f (with mean entropy %f)")
                                % alpha % m_mean(alpha) << endl;
                        m_rdirichlet.push_back(rdirichlet_t(k, alpha));
                        m_ddirichlet.push_back(ddirichlet_t(k, alpha));
                        // compute new alpha
                        alpha -= n*m_sigma(alpha)/m_dmean(alpha);
                        // increase number of components
                        m_size += 1;
                }
                // go to higher alpha values
                for (double alpha = 1.0 + m_sigma(1.0)/m_dmean(1.0); alpha < alpha_max;) {
                        cout << boost::format("adding distribution at alpha = %f (with mean entropy %f)")
                                % alpha % m_mean(alpha) << endl;
                        m_rdirichlet.push_back(rdirichlet_t(k, alpha));
                        m_ddirichlet.push_back(ddirichlet_t(k, alpha));
                        // compute new alpha
                        alpha += n*m_sigma(alpha)/m_dmean(alpha);
                        // increase number of components
                        m_size += 1;
                }
        }

        template <class Engine>
        std::vector<probability_t> operator()(Engine& eng) {
                boost::random::uniform_int_distribution<> rint(0, m_size-1);

                return m_rdirichlet[rint(eng)](eng);
        }

        double pdf(std::vector<probability_t> x) {
                double result = 0.0;

                for (size_t i = 0; i < m_size; i++) {
                        result += boost::math::pdf(m_ddirichlet[i], x);
                }
                return result/m_size;
        }
};

prob_histogram_t
approximate_distribution(size_t k, size_t minimum_counts, size_t bins)
{
        boost::random::mt19937 gen; seed_rng(gen);
        proposal_distribution_t proposal_distribution(k);
        prob_histogram_t histogram(0.0, log(k), bins);

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

        for (size_t k = 50; k <= 50; k++) {
                cerr << boost::format("Sampling entropies for cardinality %d...") % k
                     << endl;

                const prob_histogram_t histogram = approximate_distribution(k, minimum_samples, bins);

                save_table    (histogram, k);
                save_histogram(histogram, k);
        }

}
