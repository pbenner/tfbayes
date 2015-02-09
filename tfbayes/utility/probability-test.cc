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
#include <vector>
#include <iostream>

#include <boost/format.hpp>

#include <tfbayes/utility/probability.hh>

using namespace std;

typedef double real_t;
typedef probability_t<real_t> p_t;
typedef vector<p_t> p_vector_t;

void
test1(void)
{
        p_t a(0.3);
        p_t b(0.2);
        p_t c;

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
test2(real_t a, real_t b)
{
        p_t pa(a);
        p_t pb(b);

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
        cout << boost::format(" 0.25 %% 1.0: %f (%f)") % (p_t( 0.25) % 1.0) % 0.25 << endl;
        cout << boost::format(" 1.25 %% 1.0: %f (%f)") % (p_t( 1.25) % 1.0) % 0.25 << endl;
        cout << boost::format(" 2.25 %% 1.0: %f (%f)") % (p_t( 2.25) % 1.0) % 0.25 << endl;
        cout << boost::format("-0.25 %% 1.0: %f (%f)") % (p_t(-0.25) % 1.0) % 0.75 << endl;
        cout << boost::format("-1.25 %% 1.0: %f (%f)") % (p_t(-1.25) % 1.0) % 0.75 << endl;
        cout << boost::format("-2.25 %% 1.0: %f (%f)") % (p_t(-2.25) % 1.0) % 0.75 << endl;
        cout << endl;
}

void
test4()
{
        p_t a(-0.345);
        p_t b( 0.345);

        cout << "a >= 0: " << (a>=0) << endl;
        cout << "b >= 0: " << (b>=0) << endl;
        cout << "exp(a): " << exp(a) << endl;
        cout << "b^0.34: " << pow(b, real_t(0.34)) << endl;
        cout << endl;
}

void
test_accuracy1()
{
        double a = 1.0;
        double b = 1e17;

        cout << "(test 1) double:        (1.0 + 1.0e17) + -1.0e17 = "
             << (a+b)-b
             << endl;
}

void
test_accuracy2()
{
        std::vector<double> v = {1.0, 1.0e17, -1.0e17};

        cout << "(test 2) double:        (1.0 + 1.0e17) + -1.0e17 = "
             << msum(v)
             << endl;
}

void
test_accuracy3()
{
        p_t a = 1.0;
        p_t b = 1.0e17;

        cout << "(test 3) probability_t: (1.0 + 1.0e17) + -1.0e17 = "
             << (a+b)-b
             << endl;
}

void
test_accuracy4()
{
        p_vector_t v = {1.0, 1.0e17, -1.0e17};

        cout << "(test 4) probability_t: (1.0 + 1.0e17) + -1.0e17 = "
             << msum(v)
             << endl;
}

void
test_accuracy5()
{
        p_vector_t v(100000, from_log_scale(-800.0));

        cout << "(test 5) 100000*exp(-800) = "
             << std::log(kahan_sum(v))
             << endl;
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
        test4();

        test_accuracy1();
        test_accuracy2();
        test_accuracy3();
        test_accuracy4();
        test_accuracy5();
}

