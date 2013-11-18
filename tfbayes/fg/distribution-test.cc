/* Copyright (C) 2013 Philipp Benner
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

#include <distribution.hh>

using namespace std;

int
main()
{
        normal_distribution_t n(3.0,1.0);

        cout << "sufficient moment 1: " << n.moments()[0] << endl
             << "sufficient moment 2: " << n.moments()[1] << endl;

        normal_distribution_t m(n);

        cout << "sufficient moment 1: " << m.moments()[0] << endl
             << "sufficient moment 2: " << m.moments()[1] << endl;

        m = n;

        cout << "sufficient moment 1: " << m.moments()[0] << endl
             << "sufficient moment 2: " << m.moments()[1] << endl;

        cout << "density: " << n(std::vector<double>(1, 0.1))
             << endl;

        return 0;
}
