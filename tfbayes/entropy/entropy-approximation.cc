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

typedef histogram_t<double, probability_t> prob_histogram_t;

template <typename T>
void seed_rng(T& rng)
{
        struct timeval tv;
        gettimeofday(&tv, NULL);
        rng.seed(tv.tv_sec*tv.tv_usec);
}

using namespace std;

void
save_table(const prob_histogram_t& histogram, size_t k)
{
        ofstream ofs((boost::format("entropy-approximation-%d.csv") % k).str());

        BOOST_FOREACH(const double& v, histogram.midpoints()) {
                ofs << boost::format("%0.8f %0.8f") % v % histogram.pdf(v)
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

prob_histogram_t
approximate_distribution(size_t k, size_t minimum_counts, size_t bins)
{
        boost::random::mt19937 gen; seed_rng(gen);
        vector<double> alpha(k, 1.0);
        boost::random::dirichlet_distribution<double, probability_t> dist(alpha);
        prob_histogram_t histogram(0.0, log(k), bins);

        while (histogram[0] < minimum_counts) {
                histogram.add(entropy(dist(gen)));
        }
        return histogram;
}

int
main(void)
{
        vector<size_t> n = {0, 0, 100000, 100};

        for (size_t k = 2; k < n.size(); k++) {
                cerr << boost::format("Sampling entropies on simplices of dimension %d...") % k
                     << endl;

                const prob_histogram_t histogram = approximate_distribution(k, n[k], 100);

                save_table    (histogram, k);
                save_histogram(histogram, k);
        }

}
