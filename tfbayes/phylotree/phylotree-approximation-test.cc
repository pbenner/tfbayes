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
#include <sstream>
#include <sys/time.h>

#include <boost/unordered_map.hpp>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/phylotree/phylotree-gradient.hh>
#include <tfbayes/phylotree/phylotree-simplify.hh>
#include <tfbayes/phylotree/phylotree-expand.hh>
#include <tfbayes/phylotree/phylotree-approximation.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/phylotree/posterior.hh>
#include <tfbayes/phylotree/utility.hh>

using namespace std;

#define alphabet_size 4

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

boost::array<double, alphabet_size>
pt_expectation(const polynomial_t<alphabet_size>& poly,
               const exponent_t  <alphabet_size>& alpha)
{
        /* store the expectations in the array result */
        boost::array<double, alphabet_size> result;

        /* to increase the count by one the polynomial is multiplied
           with a term */
        polynomial_term_t<alphabet_size> term(1.0);
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

void test_line_search(
        const polynomial_t<alphabet_size>& result,
        const exponent_t<alphabet_size>& alpha)
{
        polynomial_t<alphabet_size> approximation;
        polynomial_t<alphabet_size> approximation_line;

        approximation      = dkl_approximate<alphabet_size>(result);
        approximation_line = dkl_line_search<alphabet_size>(approximation, result.normalize(), alpha, 100);

        cout << "h[Pa_,Pc_,Pg_,Pt_]:= "
             << approximation_line
             << endl;

        cout << "kldivl = "
             << dkl<alphabet_size>(approximation_line, result.normalize(), alpha)
             << ";"
             << endl;

}

int main(void)
{
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

        polynomial_t<alphabet_size> result = pt_polynomial_t<alphabet_size>(pt_root, observations);
        polynomial_t<alphabet_size> approximation;
        polynomial_t<alphabet_size> variational;

        cout << "f[Pa_,Pc_,Pg_,Pt_]:= " << result.normalize() << endl;

        approximation = dkl_approximate<alphabet_size>(result);

        cout << "g[Pa_,Pc_,Pg_,Pt_]:= " << approximation << endl;

        cout << "fxpt = "
             << print_expectation(pt_expectation(result, alpha))
             << endl;

        cout << "gxpt = "
             << print_expectation(pt_expectation(approximation, alpha))
             << endl;

        cout << "kldiv = "
             << dkl<alphabet_size>(approximation, result.normalize(), alpha)
             << ";"
             << endl;

        //test_line_search(result, alpha);

        variational = dkl_optimize(result, alpha);

        cout << "kldiv = "
             << dkl<alphabet_size>(variational, result.normalize(), alpha)
             << ";"
             << endl;

        return 0.0;
}
