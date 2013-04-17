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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <phylotree.hh>

using namespace std;

void test_tree1() {
        cout << "Test 1:" << endl;
        pt_leaf_t n9(9.0, "n9");
        pt_leaf_t n8(8.0, "n8");
        pt_leaf_t n7(7.0, "n7");
        pt_leaf_t n6(6.0, "n6");
        pt_node_t n5(5.0, &n8, &n9);
        pt_node_t n4(4.0, &n6, &n7);
        pt_leaf_t n3(3.0, "n3");
        pt_node_t n2(2.0, &n4, &n5);
        pt_root_t n1(&n2, &n3);

        cout << &n1 << endl;
        cout << newick_format(&n1) << endl;
}

int main(void) {
        test_tree1();

        return 0.0;
}
