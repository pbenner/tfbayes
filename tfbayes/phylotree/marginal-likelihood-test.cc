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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <tfbayes/phylotree/phylotree-expand.hh>
#include <tfbayes/phylotree/phylotree-simplify.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/phylotree/utility.hh>

using namespace std;

#define alphabet_size 5
typedef char code_t;

void test_tree1(const exponent_t<alphabet_size, code_t> alpha) {
        cout << "Test 1:" << endl;
        pt_leaf_t* n2 = new pt_leaf_t(1.0, "n2");
        pt_leaf_t* n3 = new pt_leaf_t(2.0, "n3");
        pt_root_t n1(n2, n3);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n2")->id] = 1;
        observations[n1("n3")->id] = 2;

        incomplete_expression_t incomplete_expression = pt_simplify(n1);
        polynomial_t<alphabet_size, code_t> polynomial = pt_expand<alphabet_size, code_t>(incomplete_expression, observations);
        //double result = pt_marginal_likelihood<alphabet_size, code_t>(n1, alpha);
        double result = pt_marginal_likelihood<alphabet_size, code_t>(polynomial, alpha);
        cout << result << endl
             << "-3.0468 (correct value)"
             << endl << endl;
}

void test_tree2(const exponent_t<alphabet_size, code_t> alpha) {
        cout << "Test 2:" << endl;
        pt_leaf_t* n2 = new pt_leaf_t(1.0, "n2");
        pt_leaf_t* n3 = new pt_leaf_t(2.0, "n3");
        pt_root_t n1(n2, n3);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n2")->id] = 1;
        observations[n1("n3")->id] = 1;

        double result = pt_marginal_likelihood<alphabet_size, code_t>(n1, observations, alpha);
        cout << result << endl
             << "-2.23056 (correct value)"
             << endl << endl;
}

void test_tree3(const exponent_t<alphabet_size, code_t> alpha) {
        cout << "Test 3:" << endl;
        pt_leaf_t* n5 = new pt_leaf_t(2.0, "n5");
        pt_leaf_t* n4 = new pt_leaf_t(1.0, "n4");
        pt_leaf_t* n3 = new pt_leaf_t(1.0, "n3");
        pt_node_t* n2 = new pt_node_t(0.5, n4, n5);
        pt_root_t n1(n2, n3);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;
        observations[n1("n3")->id] = 2;

        double result = pt_marginal_likelihood<alphabet_size, code_t>(n1, observations, alpha);
        cout << result << endl
             << "-4.11845 (correct value)"
             << endl << endl;
}

void test_tree4(const exponent_t<alphabet_size, code_t> alpha) {
        cout << "Test 4:" << endl;
        pt_leaf_t* n7 = new pt_leaf_t(2.0, "n7");
        pt_leaf_t* n6 = new pt_leaf_t(1.0, "n6");
        pt_leaf_t* n5 = new pt_leaf_t(2.0, "n5");
        pt_leaf_t* n4 = new pt_leaf_t(1.0, "n4");
        pt_node_t* n3 = new pt_node_t(0.5, n6, n7);
        pt_node_t* n2 = new pt_node_t(0.5, n4, n5);
        pt_root_t n1(n2, n3);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 1;
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;

        double result = pt_marginal_likelihood<alphabet_size, code_t>(n1, observations, alpha);
        cout << result << endl
             << "-3.41318 (correct value)"
             << endl << endl;
}

void test_tree5(const exponent_t<alphabet_size, code_t> alpha) {
        cout << "Test 5:" << endl;
        pt_leaf_t* n7 = new pt_leaf_t(2.0, "n7");
        pt_leaf_t* n6 = new pt_leaf_t(1.0, "n6");
        pt_leaf_t* n5 = new pt_leaf_t(2.0, "n5");
        pt_leaf_t* n4 = new pt_leaf_t(1.0, "n4");
        pt_node_t* n3 = new pt_node_t(0.5, n6, n7);
        pt_node_t* n2 = new pt_node_t(0.5, n4, n5);
        pt_root_t n1(n2, n3);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 1;
        observations[n1("n5")->id] = 2;
        observations[n1("n4")->id] = 3;

        double result = pt_marginal_likelihood<alphabet_size, code_t>(n1, observations, alpha);
        cout << result << endl
             << "-6.05774 (correct value)"
             << endl << endl;
}

void test_tree6(const exponent_t<alphabet_size, code_t> alpha) {
        cout << "Test 6:" << endl;
        pt_leaf_t* n9 = new pt_leaf_t(9.0, "n9");
        pt_leaf_t* n8 = new pt_leaf_t(8.0, "n8");
        pt_leaf_t* n7 = new pt_leaf_t(7.0, "n7");
        pt_leaf_t* n6 = new pt_leaf_t(6.0, "n6");
        pt_node_t* n5 = new pt_node_t(5.0, n8, n9);
        pt_node_t* n4 = new pt_node_t(4.0, n6, n7);
        pt_leaf_t* n3 = new pt_leaf_t(3.0, "n3");
        pt_node_t* n2 = new pt_node_t(2.0, n4, n5);
        pt_root_t n1(n2, n3);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n9")->id] = 3;
        observations[n1("n8")->id] = 2;
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 0;
        observations[n1("n3")->id] = 0;

        double result = pt_marginal_likelihood<alphabet_size, code_t>(n1, observations, alpha);
        cout << result << endl
             << "-8.1197 (correct value)"
             << endl << endl;
}

int main(void) {
        exponent_t<alphabet_size, code_t> alpha;
        alpha[0] = 1.0;
        alpha[1] = 1.0;
        alpha[2] = 1.0;
        alpha[3] = 1.0;
        alpha[4] = 0.0;

        test_tree1(alpha);
        test_tree2(alpha);
        test_tree3(alpha);
        test_tree4(alpha);
        test_tree5(alpha);
        test_tree6(alpha);

        return 0.0;
}
