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
#include <cassert>

#include <sys/time.h>

#include <tfbayes/phylotree/polynomial.hh>
#include <tfbayes/phylotree/simple-polynomial.hh>
#include <tfbayes/phylotree/utility.hh>
#include <tfbayes/uipac/alphabet.hh>

using namespace std;

#define alphabet_size 4

pt_root_t
create_tree(size_t n)
{
        assert(n >= 1);

        pt_node_t* pt_last = new pt_leaf_t(1.0);

        for (size_t i = 1; i < n; i++) {
                pt_last = new pt_node_t(1.0, pt_last, new pt_leaf_t(1.0));
        }
        return pt_root_t(pt_last, new pt_leaf_t(1.0));
}

void random_init(std::vector<alphabet_code_t>& observations)
{
        for (size_t i = 0; i < observations.size(); i++) {
                observations[i] = rand()%alphabet_size;
        }
}

double expected_terms(size_t n, bool felsenstein = false)
{
        const pt_root_t tree = create_tree(n);
        std::vector<alphabet_code_t> observations(tree.n_leaves);
        double result = 0;

        for (size_t i = 0; i < 1000; i++) {
                /* generate a new random sample */
                random_init(observations);
                /* compute the polynomial for this sample */
                if (felsenstein) {
                        pt_simple_polynomial_t<alphabet_size> poly(tree, observations);
                        /* record size of the polynomial */
                        result = (i*(double)result + (double)poly.size())/((double)i+1.0);
                }
                else {
                        polynomial_t<alphabet_size> poly = pt_likelihood<alphabet_size, alphabet_code_t, double>(
                                tree, observations);
                        /* record size of the polynomial */
                        result = (i*(double)result + (double)poly.size())/((double)i+1.0);
                }
        }

        return result;
}

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

int main(void)
{
        init();

        size_t n1 = 5;
        size_t n2 = 20;
        size_t n  = max(n1,n2);

        for (size_t i = 1; i <= n; i++) {
                /* header */
                cout << 2*i-1 << " ";
                /* felsenstein */
                if (i <= n1) {
                        cout << expected_terms(i, true) << " ";
                }
                else {
                        cout << "NA ";
                }
                /* PY */
                if (i <= n2) {
                        cout << expected_terms(i, false) << endl;
                }
                else {
                        cout << "NA" << endl;
                }
        }

        return 0.0;
}
