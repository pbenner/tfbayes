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
#include <boost/date_time/posix_time/posix_time.hpp>

#include <tfbayes/utility/histogram.hh>
#include <tfbayes/utility/probability.hh>
#include <tfbayes/utility/distribution.hh>
#include <tfbayes/entropy/entropy.hh>

template <typename T>
void seed_rng(T& rng)
{
        struct timeval tv;
        gettimeofday(&tv, NULL);
        rng.seed(tv.tv_sec*tv.tv_usec);
}

histogram_t<double, probability_t>
approximate_distribution(size_t k, size_t n, size_t bins)
{
        boost::random::mt19937 gen; seed_rng(gen);
        std::vector<double> alpha(k, 1.0);
        boost::random::dirichlet_distribution<double, probability_t> dist(alpha);
        histogram_t<double, probability_t> histogram(0.0, std::log(k), bins);

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

                const histogram_t<double, probability_t> histogram = approximate_distribution(k, n[k], 100);

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
