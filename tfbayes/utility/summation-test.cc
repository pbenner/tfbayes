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

#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>

#include <tfbayes/utility/summation.hh>

using namespace std;

typedef double real_t;
typedef vector<real_t> real_vector_t;

void
test1(void)
{
        real_vector_t v = {0.019392971714356, 0.020670102070272, 0.035717114320037, 0.061144544010206, 0.016456475532401, 0.028513948339916, 0.074321236538977, 0.030569835271186, 0.031494784585127, 0.051874277404666, 0.064300091814859, 0.050863988801977, 0.042513513743020, 0.053316344124521, 0.057613307669230, 0.056689054414379, 0.063970201966049, 0.089497305694715, 0.058533828120769, 0.092547073863335};

        cout << " sum: " << std::accumulate(v.begin(), v.end(), real_t(0.0)) << endl
             << "msum: " << msum(v)
             << endl;
}

void
test2(void)
{
        real_vector_t v = {1.0, 1.0e17, -1.0e17};

        cout << " sum: " << std::accumulate(v.begin(), v.end(), real_t(0.0)) << endl
             << "msum: " << msum(v)
             << endl;
}

int
main(void)
{
        test1();
        test2();
}
