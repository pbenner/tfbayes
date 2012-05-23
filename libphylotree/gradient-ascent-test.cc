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

#define alphabet_size 5
typedef short code_t;

using namespace std;

ostream& operator<< (ostream& o, const exponent_t<code_t, alphabet_size>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<code_t, alphabet_size>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<code_t, alphabet_size>& polynomial) {
        for (polynomial_t<code_t, alphabet_size>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + " << *it;
                }
                else {
                        o << *it;
                }
        }

        return o;
}

size_t hash_value(const exponent_t<code_t, alphabet_size>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] << 0;
        seed += (size_t)exponent[1] << 2;
        seed += (size_t)exponent[2] << 4;
        seed += (size_t)exponent[3] << 6;
        seed += (size_t)exponent[4] << 8;

        return seed;
}

void test_tree4() {
        yyparse();
        pt_root_t* pt_root = (pt_root_t*)pt_parsetree->convert();

        pt_root->init(alphabet_size);

        alignment_t<code_t> alignment("test.fa", pt_root);
//        alignment_t<code_t> alignment("test2.fa", pt_root);
        exponent_t<code_t, alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;
        alpha[4] = 1;

        cout << pt_root
             << endl;

        pt_gradient_ascent_t<code_t, alphabet_size> pt_gradient_ascent(alignment, alpha, 2.0, 0.1, 0.001);

        pt_gradient_ascent.run(1, 0.0005);

        stringstream ss;
        pt_root->print(ss, true);
        cout << ss.str()
             << endl;
}

int main(void) {
        test_tree4();

        return 0.0;
}
