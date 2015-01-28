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

#include <tfbayes/utility/probability.hh>

using namespace std;

int
main(void)
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
}
