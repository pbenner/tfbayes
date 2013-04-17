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

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

int main(void) {

        MET_INIT;
        init();

        // parse tree
        list<pt_root_t*> tree_list = parse_tree_list();
        assert(tree_list.size() == 1);
        pt_root_t* pt_root = tree_list.front();

        // random observations
        vector<code_t> observations(pt_root->n_leafs, 0);
        for (pt_node_t::id_t i = 0; i < pt_root->n_leafs; i++) {
                observations[i] = rand()%alphabet_size;
        }

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

        cout << "Parsetree:" << endl
             << pt_parsetree << endl;

        cout << "Phylogenetic tree:" << endl
             << pt_root              << endl;

        MET("Simplifying",
            incomplete_expression_t incomplete_expression = pt_simplify(pt_root));

        cout << "Simplified polynomial:" << endl
             << incomplete_expression    << endl;

        MET("Expanding",
            polynomial_t<code_t, alphabet_size> result1 = pt_expand<code_t, alphabet_size>(incomplete_expression, observations));

        cout << "Expanded polynomial:" << endl
             << result1                << endl;

        MET("Direct computation",
            polynomial_t<code_t, alphabet_size> result2 = pt_polynomial<code_t, alphabet_size>(pt_root, observations));

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
