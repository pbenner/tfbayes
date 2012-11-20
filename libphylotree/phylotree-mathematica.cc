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

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <boost/unordered_map.hpp>

#include <sys/time.h>

#include <phylotree.hh>
#include <phylotree-parser.hh>
#include <phylotree-polynomial.hh>
#include <phylotree-gradient.hh>
#include <phylotree-simplify.hh>
#include <phylotree-expand.hh>
#include <phylotree-approximation.hh>
#include <marginal-likelihood.hh>
#include <posterior.hh>
#include <utility.hh>

#include <tfbayes/exception.h>

using namespace std;

#define alphabet_size 4
typedef float code_t;

void pt_init_leaf(const boost::unordered_map<string, code_t>& observations, pt_node_t* node) {
        if (observations.find(node->name) == observations.end()) {
                node->x = rand()%alphabet_size;
        }
        else {
                node->x = (*observations.find(node->name)).second;
        }
}

void pt_init(const boost::unordered_map<string, code_t>& observations, pt_node_t* node) {

        if (node->leaf()) {
                pt_init_leaf(observations, node);
        }
        else {
                pt_init(observations, node->left);
                pt_init(observations, node->right);
        }
}

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

boost::array<double, alphabet_size>
pt_expectation(const pt_polynomial_t<code_t, alphabet_size>& poly, const exponent_t<code_t, alphabet_size>& alpha)
{
        /* store the expectations in the array result */
        boost::array<double, alphabet_size> result;

        /* to increase the count by one the polynomial is multiplied
           with a term */
        polynomial_term_t<code_t, alphabet_size> term(1.0);
        double ml = pt_marginal_likelihood(poly, alpha);

        for (size_t i = 0; i < alphabet_size; i++) {
                term.exponent()[i] = 1;
                result[i] = exp(pt_marginal_likelihood(term*poly, alpha) - ml);
                term.exponent()[i] = 0;
        }

        return result;
}

string print_expectation(const boost::array<double, alphabet_size>& expectation)
{
        stringstream ss;

        ss << "{";
        for (size_t i = 0; i < alphabet_size; i++) {
                ss << expectation[i];
                if (i != alphabet_size-1) {
                        ss << ", ";
                }
        }
        ss << "};";

        return ss.str();
}

void test(void)
{
        exponent_t<code_t, 2> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        exponent_t<double, 2> theta;
        theta[0] = 0.5;
        theta[1] = 0.5;

        polynomial_term_t<code_t, 2> term1(0.2);
        term1.exponent()[0] = 14;
        term1.exponent()[1] =  3;

        polynomial_term_t<code_t, 2> term2(0.8);
        term2.exponent()[0] = 10;
        term2.exponent()[1] =  4;

        polynomial_t<code_t, 2> poly;
        polynomial_t<code_t, 2> variational;

        poly += term1;
        poly += term2;

        cout << "RESULT: "
             << exp(poly.log_eval(theta))
             << endl
             << "ML: "
             << pt_marginal_likelihood(poly, alpha)
             << endl;

        variational = pt_approximate<code_t, 2>(poly);

        cout << kl_divergence<code_t, 2>(variational, poly, alpha)
             << endl;

}

int main(void) {

        init();
        yyparse();

        boost::unordered_map<string, code_t> observations;
        exponent_t<code_t, alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;

        pt_root_t* pt_root = (pt_root_t*)pt_parsetree->convert();

        pt_init(observations, pt_root);

        pt_polynomial_t<code_t, alphabet_size> result(pt_root);
        pt_polynomial_t<code_t, alphabet_size> variational;

        cout << "f[Pa_,Pc_,Pg_,Pt_]:= " << result.normalize() << endl;

        variational = pt_approximate<code_t, alphabet_size>(result);

        cout << "g[Pa_,Pc_,Pg_,Pt_]:= " << variational << endl;

        cout << "fxpt = " 
             << print_expectation(pt_expectation(result, alpha))
             << endl;

        cout << "gxpt = " 
             << print_expectation(pt_expectation(variational, alpha))
             << endl;

        cout << "kldiv = "
             << kl_divergence<code_t, alphabet_size>(variational, result.normalize(), alpha)
             << ";"
             << endl;

        pt_parsetree->destroy();
        pt_root->destroy();

        return 0.0;
}
