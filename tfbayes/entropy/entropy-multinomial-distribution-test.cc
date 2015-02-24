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

#include <entropy-multinomial-distribution.hh>

using namespace std;

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

        std::vector<double> counts({3.0, 1.0, 2.0, 3.0, 3.0});
        std::vector<double> theta ({0.1, 0.1, 0.1, 0.2, 0.5});

        entropy_multinomial_distribution_t<> ecat(counts.size(), 10, 20);

        cout << "join: " << pdf(ecat, theta, counts)
             << endl;
        cout << "marginal: " << marginalize(ecat, counts, gen)
             << endl;
}
