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

#include <treespace.hh>

using namespace std;

int
main(void)
{
        set<size_t> s1;
        s1.insert(3);
        s1.insert(4);
        set<size_t> s2;
        s2.insert(1);
        s2.insert(2);
        set<size_t> s3;
        s3.insert(0);
        s3.insert(5);
        set<size_t> s4;
        s4.insert(2);
        s4.insert(3);
        set<size_t> s5;
        s5.insert(4);
        s5.insert(5);
        set<size_t> s6;
        s6.insert(0);
        s6.insert(1);
        nsplit_t e1(6, s1);
        nsplit_t e2(6, s2);
        nsplit_t e3(6, s3);
        nsplit_t e4(6, s4);
        nsplit_t e5(6, s5);
        nsplit_t e6(6, s6);
        cout << e4 << endl
             << e5 << endl;
        cout << "compatible: "
             << compatible(e4, e5)
             << endl;

        return 0;
}
