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

#include <likelihood.hh>
#include <utility.hh>

using namespace std;

void test_tree1() {
        cout << "Test 1:" << endl;
        pt_leaf_t<code_t, alphabet_size> n2(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n3(2, 2.0);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        polynomial_t<code_t, alphabet_size> result = pt_likelihood(&n1);
        cout << result << endl
             << "0.950213 Pa^0 Pc^1 Pg^1 Pt^0 + 0 (correct polynomial)"
             << endl << endl;
}

void test_tree2() {
        cout << "Test 2:" << endl;
        pt_leaf_t<code_t, alphabet_size> n2(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n3(1, 2.0);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        polynomial_t<code_t, alphabet_size> result = pt_likelihood(&n1);
        cout << result << endl
             << "0.950213 Pa^0 Pc^2 Pg^0 Pt^0 + 0.0497871 Pa^0 Pc^1 Pg^0 Pt^0 + 0 (correct polynomial)"
             << endl << endl;
}

void test_tree3() {
        cout << "Test 3:" << endl;
        pt_leaf_t<code_t, alphabet_size> n5(1, 2.0);
        pt_leaf_t<code_t, alphabet_size> n4(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n3(2, 1.0);
        pt_node_t<code_t, alphabet_size> n2(-1, 0.5, &n4, &n5);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        polynomial_t<code_t, alphabet_size> result = pt_likelihood(&n1);
        cout << result << endl
             << "0.0386781 Pa^0 Pc^1 Pg^1 Pt^0 + 0.860149 Pa^0 Pc^2 Pg^1 Pt^0 + 0 (correct polynomial)"
             << endl << endl;
}

void test_tree4() {
        cout << "Test 4:" << endl;
        pt_leaf_t<code_t, alphabet_size> n7(1, 2.0);
        pt_leaf_t<code_t, alphabet_size> n6(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n5(1, 2.0);
        pt_leaf_t<code_t, alphabet_size> n4(1, 1.0);
        pt_node_t<code_t, alphabet_size> n3(-1, 0.5, &n6, &n7);
        pt_node_t<code_t, alphabet_size> n2(-1, 0.5, &n4, &n5);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        polynomial_t<code_t, alphabet_size> result = pt_likelihood(&n1);
        cout << result << endl
             << "0.0163527 Pa^0 Pc^2 Pg^0 Pt^0 + 0.000911882 Pa^0 Pc^1 Pg^0 Pt^0 + 0.842968 Pa^0 Pc^4 Pg^0 Pt^0 + 0.139768 Pa^0 Pc^3 Pg^0 Pt^0 + 0 (correct polynomial)"
             << endl << endl;
}

void test_tree5() {
        cout << "Test 5:" << endl;
        pt_leaf_t<code_t, alphabet_size> n7(1, 2.0);
        pt_leaf_t<code_t, alphabet_size> n6(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n5(2, 2.0);
        pt_leaf_t<code_t, alphabet_size> n4(3, 1.0);
        pt_node_t<code_t, alphabet_size> n3(-1, 0.5, &n6, &n7);
        pt_node_t<code_t, alphabet_size> n2(-1, 0.5, &n4, &n5);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        polynomial_t<code_t, alphabet_size> result = pt_likelihood(&n1);
        cout << result << endl
             << "0.0399154 Pa^0 Pc^1 Pg^1 Pt^1 + 0.842968 Pa^0 Pc^2 Pg^1 Pt^1 + 0 (correct polynomial)"
             << endl << endl;
}

int main(void) {
        test_tree1();
        test_tree2();
        test_tree3();
        test_tree4();
        test_tree5();

        return 0.0;
}
