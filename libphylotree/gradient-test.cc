/* Copyright (C) 2012 Philipp Benner
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <phylotree-gradient.hh>
#include <utility.hh>

using namespace std;

void test_tree1() {
        cout << "Test 1:" << endl;

        code_t observations[][2] = {
                {0, 0}, {0, 0}, {1, 0}
        };
        boost::array<double, alphabet_size> p;
        p[0] = 0.25;
        p[1] = 0.25;
        p[2] = 0.25;
        p[3] = 0.25;

        double gradient = 0;
        for (size_t i = 0; i < 3; i++) {
                pt_leaf_t n2( observations[i][0], 0.30, "n2");
                pt_leaf_t n3( observations[i][1], 0.30, "n3");
                pt_root_t n1(-1, &n2, &n3);

                pt_gradient_t<code_t, alphabet_size> result(&n1);
                gradient += result.eval(&n2, p);
        }
        cout << "gradient for n2: " << gradient
             << endl;
}

int main(void) {
        test_tree1();

        return 0.0;
}
