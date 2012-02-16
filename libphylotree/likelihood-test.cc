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

#include <phylotree-polynomial.hh>
#include <utility.hh>

using namespace std;

void test_tree1() {
        cout << "Test 1:" << endl;
        pt_leaf_t n2(1, 1.0);
        pt_leaf_t n3(2, 2.0);
        pt_root_t n1(-1, &n2, &n3);

        pt_polynomial_t<code_t, alphabet_size> result(&n1);
        cout << result << endl
             << "0.950213 Pc^1 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree2() {
        cout << "Test 2:" << endl;
        pt_leaf_t n2(1, 1.0);
        pt_leaf_t n3(1, 2.0);
        pt_root_t n1(-1, &n2, &n3);

        pt_polynomial_t<code_t, alphabet_size> result(&n1);
        cout << result << endl
             << "0.950213 Pc^2 + 0.0497871 Pc^1 (correct polynomial)"
             << endl << endl;
}

void test_tree3() {
        cout << "Test 3:" << endl;
        pt_leaf_t n5(1, 2.0);
        pt_leaf_t n4(1, 1.0);
        pt_leaf_t n3(2, 1.0);
        pt_node_t n2(-1, 0.5, &n4, &n5);
        pt_root_t n1(-1, &n2, &n3);

        pt_polynomial_t<code_t, alphabet_size> result(&n1);
        cout << result << endl
             << "0.0386781 Pc^1 Pg^1 + 0.860149 Pc^2 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree4() {
        cout << "Test 4:" << endl;
        pt_leaf_t n7(1, 2.0);
        pt_leaf_t n6(1, 1.0);
        pt_leaf_t n5(1, 2.0);
        pt_leaf_t n4(1, 1.0);
        pt_node_t n3(-1, 0.5, &n6, &n7);
        pt_node_t n2(-1, 0.5, &n4, &n5);
        pt_root_t n1(-1, &n2, &n3);

        pt_polynomial_t<code_t, alphabet_size> result(&n1);
        cout << result << endl
             << "0.0163527 Pc^2 + 0.000911882 Pc^1 + 0.842968 Pc^4 + 0.139768 Pc^3 (correct polynomial)"
             << endl << endl;
}

void test_tree5() {
        cout << "Test 5:" << endl;
        pt_leaf_t n7(1, 2.0);
        pt_leaf_t n6(1, 1.0);
        pt_leaf_t n5(2, 2.0);
        pt_leaf_t n4(3, 1.0);
        pt_node_t n3(-1, 0.5, &n6, &n7);
        pt_node_t n2(-1, 0.5, &n4, &n5);
        pt_root_t n1(-1, &n2, &n3);

        pt_polynomial_t<code_t, alphabet_size> result(&n1);
        cout << result << endl
             << "0.0399154 Pc^1 Pg^1 Pt^1 + 0.842968 Pc^2 Pg^1 Pt^1 (correct polynomial)"
             << endl << endl;
}

void test_tree6() {
        cout << "Test 6:" << endl;
        pt_leaf_t n9( 3, 9.0);
        pt_leaf_t n8( 2, 8.0);
        pt_leaf_t n7( 1, 7.0);
        pt_leaf_t n6( 0, 6.0);
        pt_node_t n5(-1, 5.0, &n8, &n9);
        pt_node_t n4(-1, 4.0, &n6, &n7);
        pt_leaf_t n3( 0, 3.0);
        pt_node_t n2(-1, 2.0, &n4, &n5);
        pt_root_t n1(-1, &n2, &n3);

        pt_polynomial_t<code_t, alphabet_size> result(&n1);
        cout << result << endl
             << "0.999997 Pa^2 Pc^1 Pg^1 Pt^1 + 3.05622e-07 Pa^1 Pc^1 Pg^1 Pt^1 (correct polynomial)"
             << endl << endl;
}

int main(void) {
        test_tree1();
        test_tree2();
        test_tree3();
        test_tree4();
        test_tree5();
        test_tree6();

        return 0.0;
}
