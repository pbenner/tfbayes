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

#include <phylotree-expand.hh>
#include <phylotree-simplify.hh>
#include <marginal-likelihood.hh>
#include <utility.hh>

using namespace std;

void test_tree1(const exponent_t<code_t, alphabet_size> alpha) {
        cout << "Test 1:" << endl;
        pt_leaf_t n2(1, 1.0);
        pt_leaf_t n3(2, 2.0);
        pt_root_t n1(-1, &n2, &n3);

        incomplete_polynomial_t incomplete_polynomial = pt_simplify(&n1);
        polynomial_t<code_t, alphabet_size> polynomial = pt_expand<code_t, alphabet_size>(incomplete_polynomial);
        //double result = pt_marginal_likelihood<code_t, alphabet_size>(&n1, alpha);
        double result = pt_marginal_likelihood<code_t, alphabet_size>(polynomial, alpha);
        cout << result << endl
             << "-3.0468 (correct value)"
             << endl << endl;
}

void test_tree2(const exponent_t<code_t, alphabet_size> alpha) {
        cout << "Test 2:" << endl;
        pt_leaf_t n2(1, 1.0);
        pt_leaf_t n3(1, 2.0);
        pt_root_t n1(-1, &n2, &n3);

        double result = pt_marginal_likelihood<code_t, alphabet_size>(&n1, alpha);
        cout << result << endl
             << "-2.23056 (correct value)"
             << endl << endl;
}

void test_tree3(const exponent_t<code_t, alphabet_size> alpha) {
        cout << "Test 3:" << endl;
        pt_leaf_t n5(1, 2.0);
        pt_leaf_t n4(1, 1.0);
        pt_leaf_t n3(2, 1.0);
        pt_node_t n2(-1, 0.5, &n4, &n5);
        pt_root_t n1(-1, &n2, &n3);

        double result = pt_marginal_likelihood<code_t, alphabet_size>(&n1, alpha);
        cout << result << endl
             << "-4.11845 (correct value)"
             << endl << endl;
}

void test_tree4(const exponent_t<code_t, alphabet_size> alpha) {
        cout << "Test 4:" << endl;
        pt_leaf_t n7(1, 2.0);
        pt_leaf_t n6(1, 1.0);
        pt_leaf_t n5(1, 2.0);
        pt_leaf_t n4(1, 1.0);
        pt_node_t n3(-1, 0.5, &n6, &n7);
        pt_node_t n2(-1, 0.5, &n4, &n5);
        pt_root_t n1(-1, &n2, &n3);

        double result = pt_marginal_likelihood<code_t, alphabet_size>(&n1, alpha);
        cout << result << endl
             << "-3.41318 (correct value)"
             << endl << endl;
}

void test_tree5(const exponent_t<code_t, alphabet_size> alpha) {
        cout << "Test 5:" << endl;
        pt_leaf_t n7(1, 2.0);
        pt_leaf_t n6(1, 1.0);
        pt_leaf_t n5(2, 2.0);
        pt_leaf_t n4(3, 1.0);
        pt_node_t n3(-1, 0.5, &n6, &n7);
        pt_node_t n2(-1, 0.5, &n4, &n5);
        pt_root_t n1(-1, &n2, &n3);

        double result = pt_marginal_likelihood<code_t, alphabet_size>(&n1, alpha);
        cout << result << endl
             << "-6.05774 (correct value)"
             << endl << endl;
}

void test_tree6(const exponent_t<code_t, alphabet_size> alpha) {
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

        double result = pt_marginal_likelihood<code_t, alphabet_size>(&n1, alpha);
        cout << result << endl
             << "-8.1197 (correct value)"
             << endl << endl;
}

int main(void) {
        exponent_t<code_t, alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;

        test_tree1(alpha);
        test_tree2(alpha);
        test_tree3(alpha);
        test_tree4(alpha);
        test_tree5(alpha);
        test_tree6(alpha);

        return 0.0;
}
