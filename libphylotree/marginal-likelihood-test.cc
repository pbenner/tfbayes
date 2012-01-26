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

#include <marginal-likelihood.hh>
#include <utility.hh>

using namespace std;

void test_tree1(const polynomial_term_t<code_t, alphabet_size> alpha) {
        cout << "Test 1:" << endl;
        pt_leaf_t<code_t, alphabet_size> n2(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n3(2, 2.0);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        double result = pt_marginal_likelihood(&n1, alpha);
        cout << result << endl
             << " (correct value)"
             << endl << endl;
}

void test_tree2(const polynomial_term_t<code_t, alphabet_size> alpha) {
        cout << "Test 2:" << endl;
        pt_leaf_t<code_t, alphabet_size> n2(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n3(1, 2.0);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        double result = pt_marginal_likelihood(&n1, alpha);
        cout << result << endl
             << " (correct value)"
             << endl << endl;
}

void test_tree3(const polynomial_term_t<code_t, alphabet_size> alpha) {
        cout << "Test 3:" << endl;
        pt_leaf_t<code_t, alphabet_size> n5(1, 2.0);
        pt_leaf_t<code_t, alphabet_size> n4(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n3(2, 1.0);
        pt_node_t<code_t, alphabet_size> n2(-1, 0.5, &n4, &n5);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        double result = pt_marginal_likelihood(&n1, alpha);
        cout << result << endl
             << " (correct value)"
             << endl << endl;
}

void test_tree4(const polynomial_term_t<code_t, alphabet_size> alpha) {
        cout << "Test 4:" << endl;
        pt_leaf_t<code_t, alphabet_size> n7(1, 2.0);
        pt_leaf_t<code_t, alphabet_size> n6(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n5(1, 2.0);
        pt_leaf_t<code_t, alphabet_size> n4(1, 1.0);
        pt_node_t<code_t, alphabet_size> n3(-1, 0.5, &n6, &n7);
        pt_node_t<code_t, alphabet_size> n2(-1, 0.5, &n4, &n5);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        double result = pt_marginal_likelihood(&n1, alpha);
        cout << result << endl
             << " (correct value)"
             << endl << endl;
}

void test_tree5(const polynomial_term_t<code_t, alphabet_size> alpha) {
        cout << "Test 5:" << endl;
        pt_leaf_t<code_t, alphabet_size> n7(1, 2.0);
        pt_leaf_t<code_t, alphabet_size> n6(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n5(2, 2.0);
        pt_leaf_t<code_t, alphabet_size> n4(3, 1.0);
        pt_node_t<code_t, alphabet_size> n3(-1, 0.5, &n6, &n7);
        pt_node_t<code_t, alphabet_size> n2(-1, 0.5, &n4, &n5);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        double result = pt_marginal_likelihood(&n1, alpha);
        cout << result << endl
             << " (correct value)"
             << endl << endl;
}

int main(void) {
        polynomial_term_t<code_t, alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;

        test_tree1(alpha);
        test_tree2(alpha);
        test_tree3(alpha);
        test_tree4(alpha);
        test_tree5(alpha);

        return 0.0;
}
