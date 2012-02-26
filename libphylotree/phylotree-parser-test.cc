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

#include <sys/time.h>

#include <phylotree.hh>
#include <phylotree-parser.hh>
#include <phylotree-polynomial.hh>
#include <phylotree-simplify.hh>
#include <phylotree-expand.hh>
#include <utility.hh>

using namespace std;

void pt_random_init(pt_node_t* node) {

        if (node->leaf()) {
                node->x = rand()%alphabet_size;
        }
        else {
                pt_random_init(node->left);
                pt_random_init(node->right);
        }
}

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

int main(void) {

        init();
        yyparse();

        pt_root_t* pt_root = (pt_root_t*)pt_parsetree->convert();

        pt_random_init(pt_root);

        cout << "Parsetree:" << endl
             << pt_parsetree << endl;

        cout << "Phylogenetic tree:" << endl
             << pt_root              << endl;


        incomplete_polynomial_t incomplete_polynomial = pt_simplify(pt_root);
// //        polynomial_t<code_t, alphabet_size> result = pt_expand<code_t, alphabet_size>(incomplete_polynomial);

        cout << "Simplified polynomial:" << endl
             << incomplete_polynomial    << endl;

        pt_parsetree->destroy();
        pt_root->destroy();

        return 0.0;
}
