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
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

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
        probability_t pdf(double value) const {
                size_t i = value == m_max ? m_n-1 : std::floor((value-m_min)/m_width);
                probability_t m = m_total;
                probability_t w = m_width;
                probability_t n;
                if (i == m_n-1) {
                        // no linear interpolation in this case
                        n = operator[](i);
                }
                else {
                        // interpolate the result
                        const double y1 = operator[](i);
                        const double y2 = operator[](i+1);
                        n = y1 + (y2 - y1)*(value - m_midpoints[i])/m_width;
                }
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
                o << "min max midpoint counts"
                  << std::endl;
                for (size_t i = 0; i < hist.m_n; i++) {
                        o << boost::format("%0.8f %0.8f %0.8f %d")
                                % min % max % hist.m_midpoints[i] % hist[i]
                          << std::endl;
                        min += hist.m_width;
                        max += hist.m_width;
                }
                return o;
        }
private:
        /* serialization */
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) {
                ar & static_cast<base_t&>(*this);
                ar & m_n;
                ar & m_min;
                ar & m_max;
                ar & m_width;
                ar & m_total;
                ar & m_midpoints;
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

int
main(void)
{
        // number of samples for each dimension
        double n[] = {0, 0, 10000000, 10000000, 50000000, 100000000};

        for (size_t k = 2; k <= 5; k++) {
                cerr << boost::format("Sampling entropies on simplices of dimension %d...") % k
                     << endl;

                const histogram_t histogram = approximate_distribution(k, n[k], 100);

                // save pdf to table
                {
                        std::ofstream ofs((boost::format("entropy-approximation-%d.csv") % k).str());

                        BOOST_FOREACH(const double& v, histogram.midpoints()) {
                                ofs << boost::format("%0.8f %0.8f") % v % histogram.pdf(v)
                                    << endl;
                        }
                }
                // save data to archive
                {
                        std::ofstream ofs((boost::format("entropy-approximation-%d.ar") % k).str());
                        boost::archive::text_oarchive oa(ofs);
                        // write class instance to archive
                        oa << histogram;
                }
        }

}
