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

#include <tfbayes/phylotree/phylotree-simplify.hh>
#include <tfbayes/phylotree/phylotree-expand.hh>
#include <tfbayes/phylotree/utility.hh>

using namespace std;

#define alphabet_size 4

void test_tree1() {
        cout << "Test 1:" << endl;
        pt_leaf_t* n2 = new pt_leaf_t(1.0, "n2");
        pt_leaf_t* n3 = new pt_leaf_t(2.0, "n3");
        pt_root_t n1(n2, n3);
        vector<alphabet_code_t> observations(n1.n_leaves, 0);
        observations[n1("n2")->id] = 1;
        observations[n1("n3")->id] = 2;

        cout << "(n1 (n2 C) (n3 G))" << endl;

        incomplete_expression_t incomplete_expression = pt_simplify(n1);
        polynomial_t<alphabet_size> result = pt_expand_t<alphabet_size>(incomplete_expression, observations);
        cout << incomplete_expression << endl
             << result << endl
             << "0.950213 Pc^1 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree2() {
        cout << "Test 2:" << endl;
        pt_leaf_t* n2 = new pt_leaf_t(1.0, "n2");
        pt_leaf_t* n3 = new pt_leaf_t(2.0, "n3");
        pt_root_t n1(n2, n3);
        vector<alphabet_code_t> observations(n1.n_leaves, 0);
        observations[n1("n2")->id] = 1;
        observations[n1("n3")->id] = 1;

        cout << "(n1 (n2 C) (n3 C))" << endl;

        incomplete_expression_t incomplete_expression = pt_simplify(n1);
        polynomial_t<alphabet_size> result = pt_expand_t<alphabet_size>(incomplete_expression, observations);
        cout << incomplete_expression << endl
             << result << endl
             << "0.950213 Pc^2 + 0.0497871 Pc^1 (correct polynomial)"
             << endl << endl;
}

void test_tree3() {
        cout << "Test 3:" << endl;
        pt_leaf_t* n5 = new pt_leaf_t(2.0, "n5");
        pt_leaf_t* n4 = new pt_leaf_t(1.0, "n4");
        pt_leaf_t* n3 = new pt_leaf_t(1.0, "n3");
        pt_node_t* n2 = new pt_node_t(0.5, n4, n5);
        pt_root_t n1(n2, n3);
        vector<alphabet_code_t> observations(n1.n_leaves, 0);
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;
        observations[n1("n3")->id] = 2;

        cout << "(n1 (n2 (n4 C) (n5 C))" << endl
             << "    (n3 G))"            << endl;

        incomplete_expression_t incomplete_expression = pt_simplify(n1);
        polynomial_t<alphabet_size> result = pt_expand_t<alphabet_size>(incomplete_expression, observations);
        cout << incomplete_expression << endl
             << result << endl
             << "0.0386781 Pc^1 Pg^1 + 0.860149 Pc^2 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree4() {
        cout << "Test 4:" << endl;
        pt_leaf_t* n7 = new pt_leaf_t(2.0, "n7");
        pt_leaf_t* n6 = new pt_leaf_t(1.0, "n6");
        pt_leaf_t* n5 = new pt_leaf_t(2.0, "n5");
        pt_leaf_t* n4 = new pt_leaf_t(1.0, "n4");
        pt_node_t* n3 = new pt_node_t(0.5, n6, n7);
        pt_node_t* n2 = new pt_node_t(0.5, n4, n5);
        pt_root_t n1(n2, n3);
        vector<alphabet_code_t> observations(n1.n_leaves, 0);
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 1;
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;

        cout << "(n1 (n2 (n4 C) (n5 C))"  << endl
             << "    (n3 (n6 C) (n7 C)))" << endl;

        incomplete_expression_t incomplete_expression = pt_simplify(n1);
        polynomial_t<alphabet_size> result = pt_expand_t<alphabet_size>(incomplete_expression, observations);
        cout << incomplete_expression << endl
             << result << endl
             << "0.0163527 Pc^2 + 0.000911882 Pc^1 + 0.842968 Pc^4 + 0.139768 Pc^3 (correct polynomial)"
             << endl << endl;
}

void test_tree5() {
        cout << "Test 5:" << endl;
        pt_leaf_t* n7 = new pt_leaf_t(2.0, "n7");
        pt_leaf_t* n6 = new pt_leaf_t(1.0, "n6");
        pt_leaf_t* n5 = new pt_leaf_t(2.0, "n5");
        pt_leaf_t* n4 = new pt_leaf_t(1.0, "n4");
        pt_node_t* n3 = new pt_node_t(0.5, n6, n7);
        pt_node_t* n2 = new pt_node_t(0.5, n4, n5);
        pt_root_t n1(n2, n3);
        vector<alphabet_code_t> observations(n1.n_leaves, 0);
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 1;
        observations[n1("n5")->id] = 2;
        observations[n1("n4")->id] = 3;

        cout << "(n1 (n2 (n4 T) (n5 G))"  << endl
             << "    (n3 (n6 C) (n7 C)))" << endl;

        incomplete_expression_t incomplete_expression = pt_simplify(n1);
        polynomial_t<alphabet_size> result = pt_expand_t<alphabet_size>(incomplete_expression, observations);
        cout << incomplete_expression << endl
             << result << endl
             << "0.0399154 Pc^1 Pg^1 Pt^1 + 0.842968 Pc^2 Pg^1 Pt^1 (correct polynomial)"
             << endl << endl;
}

void test_tree6() {
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
        vector<alphabet_code_t> observations(n1.n_leaves, 0);
        observations[n1("n9")->id] = 3;
        observations[n1("n8")->id] = 2;
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 0;
        observations[n1("n3")->id] = 0;

        cout << "(n1 (n2 (n4 (n6 A) (n7 C))"  << endl
             << "        (n5 (n8 G) (n9 T)))" << endl
             << "    (n3 A))"                 << endl;

        incomplete_expression_t incomplete_expression = pt_simplify(n1);
        polynomial_t<alphabet_size> result = pt_expand_t<alphabet_size>(incomplete_expression, observations);
        cout << incomplete_expression << endl
             << result << endl
             << "0.999997 Pa^2 Pc^1 Pg^1 Pt^1 + 3.05622e-07 Pa^1 Pc^1 Pg^1 Pt^1 (correct polynomial)"
             << endl << endl;

        for (incomplete_expression_t::const_iterator it = incomplete_expression.begin(); it != incomplete_expression.end(); it++) {
                cout << *it << ": "
                     << pt_expand_t<alphabet_size>(*it, observations) << endl;
        }
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
