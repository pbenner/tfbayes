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

#include <boost/random/mersenne_twister.hpp>

#include <tfbayes/utility/distribution.hh>

using namespace std;

template <class T>
ostream& operator<<(ostream& o, vector<T>& v) {
        o << setprecision(8)
          << fixed;
        for (size_t i = 0; i < v.size(); i++) {
                o << v[i] << " ";
        }
        return o;
}

int
main(void)
{
        boost::random::mt19937 gen;

        vector<double> alpha(5, 1.0);
        vector<probability_t> x(5, probability_t(1.0/5.0));

        boost::math  ::dirichlet_distribution<> dist(alpha);
        boost::random::dirichlet_distribution<probability_t> rdist(alpha);

        cout << boost::math::log_pdf(dist, x)
             << endl;

        for (size_t i = 0; i < 100; i++) {
                vector<probability_t> tmp = rdist(gen);
                cout << tmp << endl;
        }
}
