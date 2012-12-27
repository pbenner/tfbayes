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

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <boost/unordered_map.hpp>

#include <sys/time.h>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/phylotree/phylotree-gradient.hh>
#include <tfbayes/phylotree/phylotree-simplify.hh>
#include <tfbayes/phylotree/phylotree-expand.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/phylotree/posterior.hh>
#include <tfbayes/phylotree/utility.hh>
#include <tfbayes/exception/exception.h>

using namespace std;

#define alphabet_size 4
typedef short code_t;

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

int main(void) {

        MET_INIT;
        init();
        yyparse();

        boost::unordered_map<string, code_t> observations;
        // observations["hg19"]    = 1;
        // observations["panTro2"] = 2;
        // observations["gorGor1"] = 0;
        // observations["ponAbe2"] = 1;
        // observations["rheMac2"] = 0;
        // observations["papHam1"] = 1;
        // observations["calJac1"] = 3;
        // observations["tarSyr1"] = 0;

        exponent_t<code_t, alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;
        exponent_t<double, alphabet_size> p;
        p[0] = 0.3;
        p[1] = 0.25;
        p[2] = 0.4;
        p[3] = 0.05;

        pt_root_t* pt_root = (pt_root_t*)pt_parsetree->convert();

        pt_init(observations, pt_root);

        cout << "Parsetree:" << endl
             << pt_parsetree << endl;

        cout << "Phylogenetic tree:" << endl
             << pt_root              << endl;

        MET("Simplifying",
            incomplete_expression_t incomplete_expression = pt_simplify(pt_root));

        cout << "Simplified polynomial:" << endl
             << incomplete_expression    << endl;

        MET("Expanding",
            polynomial_t<code_t, alphabet_size> result1 = pt_expand<code_t, alphabet_size>(incomplete_expression));

        cout << "Expanded polynomial:" << endl
             << result1                << endl;

        MET("Direct computation",
            polynomial_t<code_t, alphabet_size> result2 = pt_polynomial<code_t, alphabet_size>(pt_root));

        cout << "Direct polynomial:" << endl
             << result2              << endl;

        cout << "Marginal result1: " << pt_marginal_likelihood<code_t, alphabet_size>(result1, alpha) << endl
             << "Marginal result2: " << pt_marginal_likelihood<code_t, alphabet_size>(result2, alpha) << endl
             << endl;

        cout << "Eval result1: " << result1.eval(p) << endl
             << "Eval result2: " << result2.eval(p) << endl
             << endl;

        double sum = 0;
        boost::array<double, alphabet_size> exp = pt_posterior_expectation<code_t, alphabet_size>(result2, alpha);
        for (code_t i = 0; i < alphabet_size; i++) {
                cout << "Expectation " << i << ": " << exp[i] << endl;
                sum += exp[i];
        }
        cout << "Sum: " << sum << endl;

        pt_parsetree->destroy();
        pt_root->destroy();

        return 0.0;
}