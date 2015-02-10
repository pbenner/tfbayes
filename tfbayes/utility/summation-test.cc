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

#include <tfbayes/utility/probability.hh>
#include <tfbayes/utility/summation.hh>

using namespace std;

typedef double real_t;
typedef probability_t<double> p_t;
typedef vector<real_t> real_vector_t;
typedef vector<p_t> p_vector_t;

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
void
test3()
{
        double a = 1.0;
        double b = 1e17;
        cout << "(test 1) double: (1.0 + 1.0e17) + -1.0e17 = "
             << (a+b)-b
             << endl;
}
void
test4()
{
        std::vector<double> v = {1.0, 1.0e17, -1.0e17};
        cout << "(test 2) double: (1.0 + 1.0e17) + -1.0e17 = "
             << msum(v)
             << endl;
}
void
test5()
{
        p_t a = 1.0;
        p_t b = 1.0e17;
        cout << "(test 3) probability_t: (1.0 + 1.0e17) + -1.0e17 = "
             << (a+b)-b
             << endl;
}
void
test6()
{
        p_vector_t v = {1.0, 1.0e17, -1.0e17};
        cout << "(test 4) probability_t: (1.0 + 1.0e17) + -1.0e17 = "
             << msum(v)
             << endl;
}
void
test7()
{
        vector<double> v = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7 };

        cout << "naive sum = "
             << setprecision(20)
             << std::accumulate(v.begin(), v.end(), double(0.0))
             << endl;
        cout << "     ksum = "
             << ksum(v)
             << endl;
        cout << "     msum = "
             << msum(v)
             << endl;
}
void
test8()
{
        p_vector_t v = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7 };

        cout << "naive sum = "
             << setprecision(20)
             << std::accumulate(v.begin(), v.end(), p_t(0.0))
             << endl;
        cout << "     ksum = "
             << ksum(v)
             << endl;
        cout << "     msum = "
             << msum(v)
             << endl;
}

int
main(void)
{
        test1(); cout << endl;
        test2(); cout << endl;
        test3(); cout << endl;
        test4(); cout << endl;
        test5(); cout << endl;
        test6(); cout << endl;
        test7(); cout << endl;
        test8(); cout << endl;
}
