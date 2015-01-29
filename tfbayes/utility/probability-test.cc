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

#include <cmath>

#include <boost/format.hpp>

#include <tfbayes/utility/probability.hh>

using namespace std;

void
test1(void)
{
        probability_t a(0.3);
        probability_t b(0.2);
        probability_t c;

        c = 0.4;

        cout << "a + b  = " << a+b << endl;
        cout << "a - b  = " << a-b << endl;
        cout << "a * b  = " << a*b << endl;
        cout << "a / b  = " << a/b << endl;
        cout << "a % b  = " << a%b << endl;
        cout << "a > b  = " << (a<b)  << endl;
        cout << "log(a) = " << log(a) << endl;
        cout << "c      = " << c << endl;
        cout << endl;
}

void
test2(double a, double b)
{
        probability_t pa(a);
        probability_t pb(b);

        cout << boost::format("a = %f, b = %f") % a % b << endl;
        cout << boost::format("a + b  = %f (%f)") % (pa+pb) % (a+b) << endl;
        cout << boost::format("a - b  = %f (%f)") % (pa-pb) % (a-b) << endl;
        cout << boost::format("a * b  = %f (%f)") % (pa*pb) % (a*b) << endl;
        cout << boost::format("a / b  = %f (%f)") % (pa/pb) % (a/b) << endl;
        cout << endl;
}

void
test3()
{
        cout << boost::format(" 0.25 %% 1.0: %f (%f)") % (probability_t( 0.25) % 1.0) % 0.25 << endl;
        cout << boost::format(" 1.25 %% 1.0: %f (%f)") % (probability_t( 1.25) % 1.0) % 0.25 << endl;
        cout << boost::format(" 2.25 %% 1.0: %f (%f)") % (probability_t( 2.25) % 1.0) % 0.25 << endl;
        cout << boost::format("-0.25 %% 1.0: %f (%f)") % (probability_t(-0.25) % 1.0) % 0.75 << endl;
        cout << boost::format("-1.25 %% 1.0: %f (%f)") % (probability_t(-1.25) % 1.0) % 0.75 << endl;
        cout << boost::format("-2.25 %% 1.0: %f (%f)") % (probability_t(-2.25) % 1.0) % 0.75 << endl;
}

int
main(void)
{
        test1();

        test2( 0.2,  0.3);
        test2(-0.2,  0.3);
        test2( 0.2, -0.3);
        test2(-0.2, -0.3);

        test2( 0.3,  0.2);
        test2(-0.3,  0.2);
        test2( 0.3, -0.2);
        test2(-0.3, -0.2);

        test3();
}

