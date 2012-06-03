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

#include <phylotree-polynomial.hh>
#include <phylotree-expand.hh>
#include <utility.hh>

using namespace std;

#define alphabet_size 4
typedef short code_t;

int main(void) {
        pt_leaf_t n10( 3, 10.0, "n10");
        pt_leaf_t n9( 3, 9.0, "n9");
        pt_leaf_t n8( 2, 8.0, "n8");
        pt_leaf_t n7( 1, 7.0, "n7");
        pt_leaf_t n6( 0, 6.0, "n6");
        pt_node_t n5(-1, 5.0, &n8, &n9);
        pt_node_t n4(-1, 4.0, &n6, &n7);
        pt_leaf_t n3( 0, 3.0, "n3");
        pt_node_t n2(-1, 2.0, &n4, &n5);
        pt_root_t n1(-1, &n2, &n3);

        nodeset_t nodeset1;
        nodeset1.insert(&n3);
        nodeset1.insert(&n6);
        nodeset1.insert(&n7);
        nodeset1.insert(&n8);
        nodeset1.insert(&n9);
        nodeset1.insert(&n10);

        nodeset_t nodeset2;
        nodeset2.insert(&n3);
        nodeset2.insert(&n6);

        nodeset_t nodeset3;
        nodeset3.insert(&n7);
        nodeset3.insert(&n8);
        nodeset3.insert(&n9);
        nodeset3.insert(&n10);

        polynomial_t<code_t, alphabet_size> result1;
        polynomial_t<code_t, alphabet_size> result2;

        /* phi(c; x_3, x_6, x_7, x_8, x_9, x_10) */
        result1 += pt_expand_rec<code_t, alphabet_size>(nodeset1.begin(), nodeset1.end(), 1);

        /* phi(c; x_3, x_6, x_7, x_8, x_9, x_10) =
         *   phi(c  ; x_3, x_6) phi(c  ; x_7, x_8, x_9, x_10) +
         *   phi(c  ; x_3, x_6) phi(nil; x_7, x_8, x_9, x_10) +
         *   phi(nil; x_3, x_6) phi(c  ; x_7, x_8, x_9, x_10)
         */
        result2 += pt_expand_rec<code_t, alphabet_size>(nodeset2.begin(), nodeset2.end(), 1)*
                   pt_expand_rec<code_t, alphabet_size>(nodeset3.begin(), nodeset3.end(), 1);
        result2 += pt_expand_rec<code_t, alphabet_size>(nodeset2.begin(), nodeset2.end(), 1)*
                   pt_expand_rec<code_t, alphabet_size>(nodeset3.begin(), nodeset3.end(), alphabet_size);
        result2 += pt_expand_rec<code_t, alphabet_size>(nodeset2.begin(), nodeset2.end(), alphabet_size)*
                   pt_expand_rec<code_t, alphabet_size>(nodeset3.begin(), nodeset3.end(), 1);

        cout << result1 << endl;
        cout << result2 << endl;

        return 0.0;
}
