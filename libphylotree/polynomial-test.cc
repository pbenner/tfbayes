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

#include <tfbayes/polynomial.hh>

using namespace std;

ostream& operator<< (ostream& o, const exponent_t<unsigned short, 4>& exponent) {
        o << " w^" << exponent[0]
          << " x^" << exponent[1]
          << " y^" << exponent[2]
          << " z^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<unsigned short, 4>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<unsigned short, 4>& polynomial) {
        for (polynomial_t<unsigned short, 4>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                o << polynomial_term_t<unsigned short, 4>(*it) << " + ";
        }
        o << polynomial.constant();

        return o;
}

size_t hash_value(const exponent_t<unsigned short, 4u>& exponent) {
        size_t seed = 0;
        boost::hash_combine(seed, exponent[0]);
        boost::hash_combine(seed, exponent[1]);
        boost::hash_combine(seed, exponent[2]);
        boost::hash_combine(seed, exponent[3]);
        return seed;
}

int main() {
        polynomial_term_t<unsigned short, 4> term1;
        term1.exponent()[0] = 2;
        term1.exponent()[1] = 4;
        term1.exponent()[2] = 1;
        term1.exponent()[3] = 3;
        polynomial_term_t<unsigned short, 4> term2;
        term2.exponent()[0] = 1;
        term2.exponent()[1] = 1;
        term2.exponent()[2] = 1;
        term2.exponent()[3] = 3;
        polynomial_t<unsigned short, 4> poly1;
        poly1 += term1;
        poly1 += term2;
        poly1 += term2;
        poly1 -= term2;
        poly1 -= term2;
        poly1 -= term2;
        polynomial_t<unsigned short, 4> poly2;
        poly2 += ((2*term1) * (3*term2));

        cout << "poly1: " << poly1
             << endl;
        cout << "poly2: " << poly2 + 3
             << endl;
        cout << "poly1*poly2: " << poly1 * (poly2 + 3)
             << endl;

        return 0;
}
