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

#include <boost/random/mersenne_twister.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <tfbayes/utility/probability.hh>
#include <tfbayes/entropy/entropy-multinomial-distribution.hh>

using namespace std;

typedef double real_t;
typedef probability_t<double> p_t;
typedef vector<p_t> pvector_t;

template <typename T>
void seed_rng(T& rng)
{
        struct timeval tv;
        gettimeofday(&tv, NULL);
        rng.seed(tv.tv_sec*tv.tv_usec);
}

int
main(void)
{
        boost::random::mt19937 gen; seed_rng(gen);

        std::vector<p_t   > theta  ({ 0.5,  0.1,  0.1,  0.2,  0.1});
        std::vector<real_t> counts({10.0, 10.0, 11.0, 10.0,  9.0});
        //std::vector<real_t> counts({46.0,  1.0,  1.0,  1.0,  1.0});

        // put mass on higher entropies
        entropy_multinomial_distribution_t<real_t, p_t> ecat1(theta.size(), 10,  1);
        entropy_multinomial_distribution_t<real_t, p_t> ecat2(theta.size(), 10, 10);

        marginal_entropy_distribution_t<real_t, p_t> mcat1(theta.size(), 10,  1, 100000, gen);
        marginal_entropy_distribution_t<real_t, p_t> mcat2(theta.size(), 10, 10, 100000, gen);

        cout << "joint: " << pdf(ecat1, theta, counts)
             << endl
             << "marginal: " << pdf(mcat1, counts)
             << endl;

        cout << "joint: " << pdf(ecat2, theta, counts)
             << endl
             << "marginal: " << pdf(mcat2, counts)
             << endl;
}
