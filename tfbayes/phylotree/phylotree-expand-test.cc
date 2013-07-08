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
#include <tfbayes/phylotree/phylotree-expand.hh>
#include <tfbayes/phylotree/utility.hh>

using namespace std;

#define alphabet_size 4
typedef short code_t;

int main(void) {
        pt_leaf_t n10(10.0, "n10");
        pt_leaf_t  n9(9.0, "n9");
        pt_leaf_t  n8(8.0, "n8");
        pt_leaf_t  n7(7.0, "n7");
        pt_leaf_t  n6(6.0, "n6");
        pt_node_t  n5(5.0, &n8, &n9);
        pt_node_t  n4(4.0, &n6, &n7);
        pt_leaf_t  n3(3.0, "n3");
        pt_node_t  n2(2.0, &n4, &n5);
        pt_root_t  n1(&n2, &n3);
        vector<code_t> observations(n1.n_leaves, 0);
        observations[n1("n10")->id] = 3;
        observations[n1("n9" )->id] = 3;
        observations[n1("n8" )->id] = 2;
        observations[n1("n7" )->id] = 1;
        observations[n1("n6" )->id] = 0;
        observations[n1("n3" )->id] = 0;

        leafset_t leafset1;
        leafset1.insert(&n3);
        leafset1.insert(&n6);
        leafset1.insert(&n7);
        leafset1.insert(&n8);
        leafset1.insert(&n9);
        leafset1.insert(&n10);

        leafset_t leafset2;
        leafset2.insert(&n3);
        leafset2.insert(&n6);

        leafset_t leafset3;
        leafset3.insert(&n7);
        leafset3.insert(&n8);
        leafset3.insert(&n9);
        leafset3.insert(&n10);

        polynomial_t<code_t, alphabet_size> result1;
        polynomial_t<code_t, alphabet_size> result2;

        /* phi(c; x_3, x_6, x_7, x_8, x_9, x_10) */
        result1 += pt_expand_rec<code_t, alphabet_size>(leafset1.begin(), leafset1.end(), observations, 1);

        /* phi(c; x_3, x_6, x_7, x_8, x_9, x_10) =
         *   phi(c  ; x_3, x_6) phi(c  ; x_7, x_8, x_9, x_10) +
         *   phi(c  ; x_3, x_6) phi(nil; x_7, x_8, x_9, x_10) +
         *   phi(nil; x_3, x_6) phi(c  ; x_7, x_8, x_9, x_10)
         */
        result2 += pt_expand_rec<code_t, alphabet_size>(leafset2.begin(), leafset2.end(), observations, 1)*
                   pt_expand_rec<code_t, alphabet_size>(leafset3.begin(), leafset3.end(), observations, 1);
        result2 += pt_expand_rec<code_t, alphabet_size>(leafset2.begin(), leafset2.end(), observations, 1)*
                   pt_expand_rec<code_t, alphabet_size>(leafset3.begin(), leafset3.end(), observations, alphabet_size);
        result2 += pt_expand_rec<code_t, alphabet_size>(leafset2.begin(), leafset2.end(), observations, alphabet_size)*
                   pt_expand_rec<code_t, alphabet_size>(leafset3.begin(), leafset3.end(), observations, 1);

        cout << result1 << endl;
        cout << result2 << endl;

        return 0.0;
}
