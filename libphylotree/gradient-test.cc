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

#include <alignment.hh>
#include <phylotree.hh>
#include <phylotree-parser.hh>
#include <phylotree-polynomial.hh>
#include <phylotree-gradient.hh>
#include <marginal-likelihood.hh>
#include <utility.hh>

using namespace std;

void test_tree1() {
        cout << "Test 1:" << endl;
        pt_leaf_t n2( 2, 0.30, "n2");
        pt_leaf_t n3( 2, 0.30, "n3");
        pt_root_t n1(-1, &n2, &n3);

        pt_symbolic_gradient_t  <code_t, alphabet_size> result1(&n1);
        pt_polynomial_t<code_t, alphabet_size> result2(&n1);

        cout << result1.normalization() << endl;

        boost::array<double, alphabet_size> p;
        p[0] = 0.25;
        p[1] = 0.25;
        p[2] = 0.25;
        p[3] = 0.25;

        cout << result1.normalization().eval(p).eval()
             << endl;
        cout << result2.eval(p)
             << endl
             << endl;
}

void test_tree2() {
        cout << "Test 2:" << endl;
        pt_leaf_t n9( 3, 9.0,           "n9");
        pt_leaf_t n8( 2, 8.0,           "n8");
        pt_leaf_t n7( 1, 7.0,           "n7");
        pt_leaf_t n6( 0, 6.0,           "n6");
        pt_node_t n5(-1, 5.0, &n8, &n9, "n5");
        pt_node_t n4(-1, 4.0, &n6, &n7, "n4");
        pt_leaf_t n3( 0, 3.0,           "n3");
        pt_node_t n2(-1, 2.0, &n4, &n5, "n2");
        pt_root_t n1(-1, &n2, &n3,      "n1");

        pt_symbolic_gradient_t  <code_t, alphabet_size> result1(&n1);
        pt_gradient_t  <code_t, alphabet_size> result2(&n1);
        pt_polynomial_t<code_t, alphabet_size> result3(&n1);

        cout << "symbolic likelihood: "
             << result1.normalization()
             << endl << endl;
        cout << "symbolic gradient: "
             << result1[&n7]
             << endl << endl;

        boost::array<double, alphabet_size> p;
        p[0] = 0.25;
        p[1] = 0.25;
        p[2] = 0.25;
        p[3] = 0.25;

        cout << "eval: "
             << result1.normalization().eval(p).eval()
             << endl;
        cout << "eval: "
             << result2.normalization().eval(p)
             << endl;
        cout << "eval: "
             << result3.eval(p)
             << endl;
}

void test_tree3() {
        cout << "Test 3:" << endl;
        pt_leaf_t n9( 3, 9.0,           "n9");
        pt_leaf_t n8( 2, 8.0,           "n8");
        pt_leaf_t n7( 1, 7.0,           "n7");
        pt_leaf_t n6( 0, 6.0,           "n6");
        pt_node_t n5(-1, 5.0, &n8, &n9, "n5");
        pt_node_t n4(-1, 4.0, &n6, &n7, "n4");
        pt_leaf_t n3( 0, 3.0,           "n3");
        pt_node_t n2(-1, 2.0, &n4, &n5, "n2");
        pt_root_t n1(-1, &n2, &n3,      "n1");

        mutation_tree_t tree1(pmut_t(&n9, false));
        mutation_tree_t tree2(pmut_t(&n8, true));

        tree1 += tree2;
        tree1 *= tree2;

        cout <<  tree1 << endl;
        cout << -tree1 << endl;
        cout <<  tree1 << endl;

        polynomial_t<code_t, alphabet_size, mutation_tree_t> poly1(tree1);
        polynomial_t<code_t, alphabet_size, mutation_tree_t> poly2(1.0);
        cout << poly1*poly1+poly1 << endl;
        cout << poly1*poly1*poly2 << endl;
}

void gradient_ascent(alignment_t<code_t>& alignment, const exponent_t<code_t, alphabet_size>& alpha) {
        pt_node_t::nodes_t nodes = alignment.tree->get_nodes();
        boost::unordered_map<pt_node_t*, double> sum;
        gamma_gradient_t gamma_gradient(0.2, 0.4);

        // initialize sum with gradients of the gamma distribution
        for (pt_node_t::nodes_t::iterator it = nodes.begin(); it != nodes.end(); it++) {
                sum[*it] = gamma_gradient.eval((*it)->d);
        }

        // loop through the alignment
        for (alignment_t<code_t>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                pt_gradient_t<code_t, alphabet_size> gradient(alignment.tree);

                double norm = 0;
                // loop over monomials
                for (polynomial_t<code_t, alphabet_size>::const_iterator ut = gradient.normalization().begin();
                     ut != gradient.normalization().end(); ut++) {
                        norm += ut->coefficient()*exp(mbeta_log(ut->exponent(), alpha));
                }
                // loop over nodes
                for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                        double result = 0;
                        for (polynomial_t<code_t, alphabet_size>::const_iterator ut = gradient[*is].begin();
                             ut != gradient[*is].end(); ut++) {
                                result += ut->coefficient()*exp(mbeta_log(ut->exponent(), alpha));
                        }
                        sum[*is] += result/norm;
                        cout << "result: "
                             << result/norm
                             << endl;
                }
        }
}

void test_tree4() {
        yyparse();
        pt_root_t* pt_root = (pt_root_t*)pt_parsetree->convert();

        pt_root->init(alphabet_size);

        alignment_t<code_t> alignment("test.fa", pt_root);

        cout << pt_root
             << endl;

        exponent_t<code_t, alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;
        gradient_ascent(alignment, alpha);
}

int main(void) {
        test_tree1();
        test_tree2();
        test_tree3();
        test_tree4();

        return 0.0;
}
