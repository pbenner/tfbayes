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

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <tfbayes/utility/probability.hh>
#include <tfbayes/utility/distribution.hh>
#include <tfbayes/entropy/entropy.hh>

class histogram_t : public std::vector<double>
{
        typedef std::vector<double> base_t;
public:
        histogram_t(double min, double max, size_t n)
                : base_t     (n, 0.0),
                  m_n        (n),
                  m_min      (min),
                  m_max      (max),
                  m_width    ((max-min)/n),
                  m_total    (0.0),
                  m_midpoints(n, 0.0) {
                for (size_t i = 0; i < n; i++) {
                        m_midpoints[i] = m_width/2.0 + i*m_width;
                }
        }

        void add(double value) {
                size_t i = value == m_max ? m_n-1 : std::floor((value-m_min)/m_width);
                assert(value <= m_max);
                assert(i < size());
                operator[](i) += 1.0;
                m_total += 1.0;
        }
        probability_t pdf(double value) {
                size_t i = value == m_max ? m_n-1 : std::floor((value-m_min)/m_width);
                probability_t n = operator[](i);
                probability_t m = m_total;
                probability_t w = m_width;
                return n/(m*w);
        }

        const size_t& n() const {
                return m_n;
        }
        const double& min() const {
                return m_min;
        }
        const double& max() const {
                return m_max;
        }
        const double& width() const {
                return m_width;
        }
        const std::vector<double>& midpoints() const {
                return m_midpoints;
        }

        friend
        std::ostream& operator<<(std::ostream& o, const histogram_t& hist) {
                double min = hist.m_min;
                double max = hist.m_min+hist.m_width;
                for (size_t i = 0; i < hist.m_n; i++) {
                        o << boost::format("[%0.5f, %0.5f): %d")
                                % min % max % hist[i]
                          << std::endl;
                        min += hist.m_width;
                        max += hist.m_width;
                }
                return o;
        }

protected:
        size_t m_n;
        double m_min;
        double m_max;
        double m_width;
        double m_total;
        std::vector<double> m_midpoints;
};

histogram_t
approximate_distribution(size_t k, size_t n, size_t bins)
{
        boost::random::mt19937 gen;
        std::vector<double> alpha(k, 1.0);
        boost::random::dirichlet_distribution<double, probability_t> dist(alpha);
        histogram_t histogram(0.0, std::log(k), bins);

        for (size_t i = 0; i < n; i++) {
                histogram.add(entropy(dist(gen)));
        }
        return histogram;
}

using namespace std;

void
test_histogram(void)
{
        histogram_t hist(0, log(3), 5);

        hist.add(0.0);
        hist.add(0.1);
        hist.add(0.2);
        hist.add(0.21);
        hist.add(0.22);
        hist.add(log(3));

        cout << hist << endl;
}

int
main(void)
{
        size_t k = 3;
        double b = 100;
        histogram_t histogram = approximate_distribution(k, 10000000, b);

        BOOST_FOREACH(const double& v, histogram.midpoints()) {
                cout << boost::format("%0.8f %0.8f") % v % histogram.pdf(v)
                     << endl;
        }
}
