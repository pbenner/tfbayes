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

#define alphabet_size 4
typedef short code_t;

void test_tree1() {
        cout << "Test 1:" << endl;
        pt_leaf_t n2(1.0, "n2");
        pt_leaf_t n3(2.0, "n3");
        pt_root_t n1(&n2, &n3);
        vector<code_t> observations(n1.n_leafs, 0);
        observations[n1("n2")->id] = 1;
        observations[n1("n3")->id] = 2;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(&n1, observations);
        cout << result << endl
             << "0.950213 Pc^1 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree2() {
        cout << "Test 2:" << endl;
        pt_leaf_t n2(1.0, "n2");
        pt_leaf_t n3(2.0, "n3");
        pt_root_t n1(&n2, &n3);
        vector<code_t> observations(n1.n_leafs, 0);
        observations[n1("n2")->id] = 1;
        observations[n1("n3")->id] = 1;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(&n1, observations);
        cout << result << endl
             << "0.950213 Pc^2 + 0.0497871 Pc^1 (correct polynomial)"
             << endl << endl;
}

void test_tree3() {
        cout << "Test 3:" << endl;
        pt_leaf_t n5(2.0, "n5");
        pt_leaf_t n4(1.0, "n4");
        pt_leaf_t n3(1.0, "n3");
        pt_node_t n2(0.5, &n4, &n5);
        pt_root_t n1(&n2, &n3);
        vector<code_t> observations(n1.n_leafs, 0);
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;
        observations[n1("n3")->id] = 2;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(&n1, observations);
        cout << result << endl
             << "0.0386781 Pc^1 Pg^1 + 0.860149 Pc^2 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree4() {
        cout << "Test 4:" << endl;
        pt_leaf_t n7(2.0, "n7");
        pt_leaf_t n6(1.0, "n6");
        pt_leaf_t n5(2.0, "n5");
        pt_leaf_t n4(1.0, "n4");
        pt_node_t n3(0.5, &n6, &n7);
        pt_node_t n2(0.5, &n4, &n5);
        pt_root_t n1(&n2, &n3);
        vector<code_t> observations(n1.n_leafs, 0);
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 1;
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(&n1, observations);
        cout << result << endl
             << "0.0163527 Pc^2 + 0.000911882 Pc^1 + 0.842968 Pc^4 + 0.139768 Pc^3 (correct polynomial)"
             << endl << endl;
}

void test_tree5() {
        cout << "Test 5:" << endl;
        pt_leaf_t n7(2.0, "n7");
        pt_leaf_t n6(1.0, "n6");
        pt_leaf_t n5(2.0, "n5");
        pt_leaf_t n4(1.0, "n4");
        pt_node_t n3(0.5, &n6, &n7);
        pt_node_t n2(0.5, &n4, &n5);
        pt_root_t n1(&n2, &n3);
        vector<code_t> observations(n1.n_leafs, 0);
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 1;
        observations[n1("n5")->id] = 2;
        observations[n1("n4")->id] = 3;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(&n1, observations);
        cout << result << endl
             << "0.0399154 Pc^1 Pg^1 Pt^1 + 0.842968 Pc^2 Pg^1 Pt^1 (correct polynomial)"
             << endl << endl;
}

void test_tree6() {
        cout << "Test 6:" << endl;
        pt_leaf_t n9(9.0, "n9");
        pt_leaf_t n8(8.0, "n8");
        pt_leaf_t n7(7.0, "n7");
        pt_leaf_t n6(6.0, "n6");
        pt_node_t n5(5.0, &n8, &n9);
        pt_node_t n4(4.0, &n6, &n7);
        pt_leaf_t n3(3.0, "n3");
        pt_node_t n2(2.0, &n4, &n5);
        pt_root_t n1(&n2, &n3);
        vector<code_t> observations(n1.n_leafs, 0);
        observations[n1("n9")->id] = 3;
        observations[n1("n8")->id] = 2;
        observations[n1("n7")->id] = 1;
        observations[n1("n6")->id] = 0;
        observations[n1("n3")->id] = 0;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(&n1, observations);
        cout << result << endl
             << "0.999997 Pa^2 Pc^1 Pg^1 Pt^1 + 3.05622e-07 Pa^1 Pc^1 Pg^1 Pt^1 (correct polynomial)"
             << endl << endl;

}

void test_tree7() {
        cout << "Test 7:" << endl;
        // first tree
        pt_leaf_t na(0.1, "A");
        pt_leaf_t nb(0.3, "B");
        pt_leaf_t nc(0.5, "C");
        pt_leaf_t nd(0.6, "D");
        pt_node_t n3(0.4, &nc, &nd);
        pt_node_t n2(0.2, &nb, &n3);
        pt_root_t n1(     &na, &n2);
        vector<code_t> observations1(n1.n_leafs, 0);
        observations1[n1("A")->id] = 0;
        observations1[n1("B")->id] = 1;
        observations1[n1("C")->id] = 0;
        observations1[n1("D")->id] = 0;
        // second tree
        pt_leaf_t ma(0.3, "A");
        pt_leaf_t mb(0.3, "B");
        pt_leaf_t mc(0.5, "C");
        pt_leaf_t md(0.6, "D");
        pt_node_t m3(0.2, &ma, &mb);
        pt_node_t m2(0.2, &mc, &md);
        pt_root_t m1(     &m2, &m3);
        vector<code_t> observations2(m1.n_leafs, 0);
        observations2[m1("A")->id] = 0;
        observations2[m1("B")->id] = 1;
        observations2[m1("C")->id] = 0;
        observations2[m1("D")->id] = 0;

        polynomial_t<code_t, alphabet_size> result1 = pt_polynomial<code_t, alphabet_size>(&n1, observations1);
        polynomial_t<code_t, alphabet_size> result2 = pt_polynomial<code_t, alphabet_size>(&m1, observations2);

        cout << result1 << endl
             << result2 << endl
             << endl;
} 

void test_tree8() {
        cout << "Test 8:" << endl;
        pt_leaf_t n5(2.0, "n5");
        pt_leaf_t n4(1.0, "n4");
        pt_leaf_t n3(1.0, "n3");
        pt_node_t n2(0.5, &n4, &n5);
        pt_leaf_t n0(3.0, "n0");
        pt_root_t n1(&n2, &n3, &n0);
        vector<code_t> observations(n1.n_leafs, 0);
        observations[n1("n5")->id] = 1;
        observations[n1("n4")->id] = 1;
        observations[n1("n3")->id] = 2;
        observations[n1("n0")->id] = 1;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(&n1, observations);
        cout << result << endl
             << "0.839524 Pc^3 Pg^1 + 0.000950355 Pc^1 Pg^1 + 0.0450738 Pc^2 Pg^1 (correct polynomial)"
             << endl << endl;
}

void test_tree9() {
        cout << "Test 9:" << endl;
        pt_leaf_t n5(2.0, "n5");
        pt_leaf_t n4(1.0, "n4");
        pt_leaf_t n3(1.0, "n3");
        pt_leaf_t n0(3.0, "n0");
        pt_node_t n1(0.5, &n0, &n3);
        pt_root_t n2(&n1, &n5, &n4);
        vector<code_t> observations(n2.n_leafs, 0);
        observations[n2("n5")->id] = 1;
        observations[n2("n4")->id] = 1;
        observations[n2("n3")->id] = 2;
        observations[n2("n0")->id] = 1;

        polynomial_t<code_t, alphabet_size> result = pt_polynomial<code_t, alphabet_size>(&n2, observations);
        cout << result << endl
             << "0.839524 Pc^3 Pg^1 + 0.000950355 Pc^1 Pg^1 + 0.0450738 Pc^2 Pg^1 (correct polynomial)"
             << endl << endl;
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

        return 0.0;
}
