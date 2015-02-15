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
#include <tfbayes/entropy/entropy-approximation-2.hh>
#include <tfbayes/entropy/entropy-approximation-3.hh>
#include <tfbayes/entropy/entropy-approximation-4.hh>
#include <tfbayes/entropy/entropy-approximation-5.hh>
#include <tfbayes/entropy/entropy.hh>

using namespace std;

typedef histogram_t<double, probability_t<> > prob_histogram_t;

int
main(void)
{
        size_t k = 3;
        prob_histogram_t histogram;

        switch (k) {
        case 2:
                histogram = prob_histogram_t(0.0, log(k), entropy_histogram_2_counts);
                break;
        case 3:
                histogram = prob_histogram_t(0.0, log(k), entropy_histogram_3_counts);
                break;
        case 4:
                histogram = prob_histogram_t(0.0, log(k), entropy_histogram_4_counts);
                break;
        case 5:
                histogram = prob_histogram_t(0.0, log(k), entropy_histogram_5_counts);
                break;
        }

        for (double v = 0.0; v <= log(k); v += 0.001) {
                cout << boost::format("%f %f") % v % pdf(histogram, v)
                     << endl;
        }
}
