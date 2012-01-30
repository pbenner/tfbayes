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

#include <phylotree.hh>
#include <utility.hh>

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
        if (polynomial.constant()) {
                o << " + " << polynomial.constant();
        }

        return o;
}

size_t hash_value(const exponent_t<code_t, alphabet_size>& exponent) {
        size_t seed = 0;
        boost::hash_combine(seed, exponent[0]);
        boost::hash_combine(seed, exponent[1]);
        boost::hash_combine(seed, exponent[2]);
        boost::hash_combine(seed, exponent[3]);
        return seed;
}
