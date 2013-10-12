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

#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/phylotree/utility.hh>

using namespace std;

#define alphabet_size 5
typedef short code_t;

void test_tree1() {
        cout << "Test 1:" << endl;
        pt_leaf_t* n2 = new pt_leaf_t(1.0, "n2");
        pt_leaf_t* n3 = new pt_leaf_t(2.0, "n3");
        pt_root_t n1(n2, n3);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n2")->id] = 1;
        observations[n1("n3")->id] = 2;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n1, observations);
        cout << result << endl
             << "0.950213 Pc^1 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree2() {
        cout << "Test 2:" << endl;
        pt_leaf_t* n2 = new pt_leaf_t(1.0, "n2");
        pt_leaf_t* n3 = new pt_leaf_t(2.0, "n3");
        pt_root_t n1(n2, n3);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n2")->id] = 1;
        observations[n1("n3")->id] = 1;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n1, observations);
        cout << result << endl
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
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;
        observations[n1("n3")->id] = 2;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n1, observations);
        cout << result << endl
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
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 1;
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n1, observations);
        cout << result << endl
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
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 1;
        observations[n1("n5")->id] = 2;
        observations[n1("n4")->id] = 3;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n1, observations);
        cout << result << endl
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
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n9")->id] = 3;
        observations[n1("n8")->id] = 2;
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 0;
        observations[n1("n3")->id] = 0;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n1, observations);
        cout << result << endl
             << "0.999997 Pa^2 Pc^1 Pg^1 Pt^1 + 3.05622e-07 Pa^1 Pc^1 Pg^1 Pt^1 (correct polynomial)"
             << endl << endl;

}

void test_tree7() {
        cout << "Test 7:" << endl;
        // first tree
        pt_leaf_t* na = new pt_leaf_t(0.1, "A");
        pt_leaf_t* nb = new pt_leaf_t(0.3, "B");
        pt_leaf_t* nc = new pt_leaf_t(0.5, "C");
        pt_leaf_t* nd = new pt_leaf_t(0.6, "D");
        pt_node_t* n3 = new pt_node_t(0.4, nc, nd);
        pt_node_t* n2 = new pt_node_t(0.2, nb, n3);
        pt_root_t n1(na, n2);
        vector<code_t> observations1(n1.n_leaves, 0);
        observations1[n1("A")->id] = 0;
        observations1[n1("B")->id] = 1;
        observations1[n1("C")->id] = 0;
        observations1[n1("D")->id] = 0;
        // second tree
        pt_leaf_t* ma = new pt_leaf_t(0.3, "A");
        pt_leaf_t* mb = new pt_leaf_t(0.3, "B");
        pt_leaf_t* mc = new pt_leaf_t(0.5, "C");
        pt_leaf_t* md = new pt_leaf_t(0.6, "D");
        pt_node_t* m3 = new pt_node_t(0.2, ma, mb);
        pt_node_t* m2 = new pt_node_t(0.2, mc, md);
        pt_root_t m1(m2, m3);
        vector<code_t> observations2(m1.n_leaves, 0);
        observations2[m1("A")->id] = 0;
        observations2[m1("B")->id] = 1;
        observations2[m1("C")->id] = 0;
        observations2[m1("D")->id] = 0;

        polynomial_t<code_t, alphabet_size> result1 = pt_polynomial<code_t, alphabet_size>(n1, observations1);
        polynomial_t<code_t, alphabet_size> result2 = pt_polynomial<code_t, alphabet_size>(m1, observations2);

        cout << result1 << " (first tree)" << endl
             << result2 << " (second tree with root at different location)" << endl
             << endl;
} 

void test_tree8() {
        cout << "Test 8:" << endl;
        pt_leaf_t* n3 = new pt_leaf_t(0.1, "n3");
        pt_leaf_t* n2 = new pt_leaf_t(0.2, "n2");
        pt_leaf_t* n0 = new pt_leaf_t(0.3, "n0");
        pt_root_t n1(n2, n3, n0);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n0")->id] = 2;
        observations[n1("n2")->id] = 1;
        observations[n1("n3")->id] = 1;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n1, observations);
        cout << result << endl
             << "0.192007 Pc^1 Pg^1 + 0.0799543 Pc^2 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree9() {
        cout << "Test 9:" << endl;
        pt_leaf_t* n5 = new pt_leaf_t(2.0, "n5");
        pt_leaf_t* n4 = new pt_leaf_t(1.0, "n4");
        pt_leaf_t* n3 = new pt_leaf_t(1.0, "n3");
        pt_node_t* n2 = new pt_node_t(0.5, n4, n5);
        pt_leaf_t* n0 = new pt_leaf_t(3.0, "n0");
        pt_root_t n1(     n2, n3, n0);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;
        observations[n1("n3")->id] = 2;
        observations[n1("n0")->id] = 1;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n1, observations);
        cout << result << endl
             << "0.839524 Pc^3 Pg^1 + 0.000950355 Pc^1 Pg^1 + 0.0450738 Pc^2 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree10() {
        cout << "Test 10:" << endl;
        pt_leaf_t* n5 = new pt_leaf_t(2.0, "n5");
        pt_leaf_t* n4 = new pt_leaf_t(1.0, "n4");
        pt_leaf_t* n3 = new pt_leaf_t(1.0, "n3");
        pt_leaf_t* n0 = new pt_leaf_t(3.0, "n0");
        pt_node_t* n1 = new pt_node_t(0.5, n0, n3);
        pt_root_t n2(n1, n5, n4);
        vector<code_t> observations(n2.n_leaves, 0);
        observations[n2("n5")->id] = 1;
        observations[n2("n4")->id] = 1;
        observations[n2("n3")->id] = 2;
        observations[n2("n0")->id] = 1;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n2, observations);
        cout << result << endl
             << "0.839524 Pc^3 Pg^1 + 0.000950355 Pc^1 Pg^1 + 0.0450738 Pc^2 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree11() {
        cout << "Test 11:" << endl;
        pt_leaf_t* n0 = new pt_leaf_t(0.7, "n0");
        pt_leaf_t* n1 = new pt_leaf_t(0.3, "n1");
        pt_leaf_t* n2 = new pt_leaf_t(0.1, "n2");
        pt_leaf_t* n3 = new pt_leaf_t(0.2, "n3");
        pt_leaf_t* n4 = new pt_leaf_t(0.6, "n4");
        pt_node_t* n5 = new pt_node_t(0.4, n2, n3);
        pt_node_t* n6 = new pt_node_t(0.5, n1, n5);
        pt_root_t n7(n6, n4, n0);
        vector<code_t> observations(n7.n_leaves, 0);
        observations[n7("n0")->id] = 1;
        observations[n7("n1")->id] = 0;
        observations[n7("n2")->id] = 1;
        observations[n7("n3")->id] = 1;
        observations[n7("n4")->id] = 3;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(n7, observations);
        cout << result << endl
             << "0.0665841 Pa^1 Pc^3 Pt^1 + 0.183048 Pa^1 Pc^2 Pt^1 + 0.0174905 Pa^1 Pc^1 Pt^1 (correct polynomial)"
             << endl << endl;
}

void test_tree12() {
        cout << "Test 12:" << endl;
        // first tree
        pt_leaf_t* n5 = new pt_leaf_t(2.0, "n5");
        pt_leaf_t* n4 = new pt_leaf_t(1.0, "n4");
        pt_leaf_t* n3 = new pt_leaf_t(1.0, "n3");
        pt_node_t* n2 = new pt_node_t(0.5, n4, n5);
        pt_root_t n1(n2, n3);
        vector<code_t> observations1(n1.n_leaves, 0);
        observations1[n1("n5")->id] = -1;
        observations1[n1("n4")->id] =  1;
        observations1[n1("n3")->id] =  2;
        // second tree
        pt_leaf_t* m2 = new pt_leaf_t(1.5, "m2");
        pt_leaf_t* m3 = new pt_leaf_t(1.0, "m3");
        pt_root_t m1(m2, m3);
        vector<code_t> observations2(m1.n_leaves, 0);
        observations2[m1("m2")->id] = 1;
        observations2[m1("m3")->id] = 2;

        polynomial_t<code_t, alphabet_size> result1 = pt_polynomial<code_t, alphabet_size>(n1, observations1);
        polynomial_t<code_t, alphabet_size> result2 = pt_polynomial<code_t, alphabet_size>(m1, observations2);
        cout << result1 << " (observation with gap '-')" << endl
             << result2 << " (tree with leaf removed)"   << endl
             << endl;
}

int main(void) {
        test_tree1();
        test_tree2();
        test_tree3();
        test_tree4();
        test_tree5();
        test_tree6();
        test_tree7();
        test_tree8();
        test_tree9();
        test_tree10();
        test_tree11();
        test_tree12();

        return 0.0;
}
