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
        cout << "Test 1:" << endl << endl;

        pt_leaf_t* n9 = new pt_leaf_t(9.0, "n9");
        pt_leaf_t* n8 = new pt_leaf_t(8.0, "n8");
        pt_leaf_t* n7 = new pt_leaf_t(7.0, "n7");
        pt_leaf_t* n6 = new pt_leaf_t(6.0, "n6");
        pt_node_t* n5 = new pt_node_t(5.0, n8, n9);
        pt_node_t* n4 = new pt_node_t(4.0, n6, n7);
        pt_leaf_t* n3 = new pt_leaf_t(3.0, "n3");
        pt_node_t* n2 = new pt_node_t(2.0, n4, n5);
        pt_leaf_t* outgroup = new pt_leaf_t(10.0, "outgroup");
        pt_root_t n1(n2, n3, outgroup);

        cout << n1 << endl;
        cout << newick_format(n1) << endl;
}

void test_tree2() {
        cout << "Test 2:" << endl << endl;

        pt_leaf_t* n5 = new pt_leaf_t(7.0, "(1)");
        pt_leaf_t* n4 = new pt_leaf_t(6.0, "(2)");
        pt_leaf_t* n3 = new pt_leaf_t(3.0, "(3)");
        pt_node_t* n2 = new pt_node_t(2.0, n4, n5);
        pt_root_t n1(n2, n3);

        cout << n1 << endl; n2->move_a();
        cout << n1 << endl; n2->move_a();
        cout << n1 << endl; n2->move_b();
        cout << n1 << endl; n2->move_b();
        cout << n1 << endl;
}

void test_tree3() {
        cout << "Test 3:" << endl << endl;

        pt_leaf_t* n5 = new pt_leaf_t(7.0, "(1)");
        pt_leaf_t* n4 = new pt_leaf_t(6.0, "(2)");
        pt_leaf_t* n3 = new pt_leaf_t(3.0, "(3)");
        pt_node_t* n2 = new pt_node_t(2.0, n4, n5);
        pt_root_t n1(n3, n2);

        cout << n1 << endl; n2->move_a();
        cout << n1 << endl; n2->move_a();
        cout << n1 << endl; n2->move_b();
        cout << n1 << endl; n2->move_b();
        cout << n1 << endl;
}

int main(void) {
        test_tree1(); cout << endl;
        test_tree2(); cout << endl;
        test_tree3(); cout << endl;

        return 0.0;
}
