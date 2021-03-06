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
#include <tfbayes/phylotree/parser.hh>
#include <tfbayes/phylotree/polynomial.hh>
#include <tfbayes/phylotree/gradient.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/phylotree/posterior.hh>
#include <tfbayes/phylotree/tree-reduction.hh>
#include <tfbayes/phylotree/utility.hh>

using namespace std;

#define alphabet_size 4

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

int main(void) {

        init();

        // parse tree
        list<pt_root_t> tree_list = parse_tree_list();
        assert(tree_list.size() == 1);
        pt_root_t& pt_root = tree_list.front();

        // random observations
        vector<alphabet_code_t> observations(pt_root.n_leaves, 0);
        for (pt_node_t::id_t i = 0; i < pt_root.n_leaves; i++) {
                observations[i] = rand() % alphabet_size;
        }

        exponent_t<alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;
        exponent_t<alphabet_size, double> p;
        p[0] = 0.3;
        p[1] = 0.25;
        p[2] = 0.4;
        p[3] = 0.05;

        cout << "Phylogenetic tree:"   << endl
             << newick_format(pt_root) << endl;

        polynomial_t<alphabet_size> result2 = pt_likelihood<alphabet_size, alphabet_code_t, double>(
                pt_root, observations);

        cout << "Direct polynomial:" << endl
             << result2              << endl;

        cout //<< "Marginal result1: " << pt_marginal_likelihood<alphabet_size>(result1, alpha) << endl
             << "Marginal result2: " << pt_marginal_likelihood(result2, alpha) << endl
             << endl;

        cout //<< "Eval result1: " << result1.eval(p) << endl
             << "Eval result2: " << result2.eval(p) << endl
             << endl;

        double sum = 0;
        boost::array<double, alphabet_size> exp = pt_posterior_expectation<alphabet_size>(
                result2, alpha);
        for (size_t i = 0; i < alphabet_size; i++) {
                cout << "Expectation " << i << ": " << exp[i] << endl;
                sum += exp[i];
        }
        cout << "Sum: " << sum << endl;

        return 0.0;
}
