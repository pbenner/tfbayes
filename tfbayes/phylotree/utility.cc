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

#include <tfbayes/phylotree/utility.hh>

using namespace std;

ostream& operator<< (ostream& o, const exponent_t<short, 4>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<short, 4>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<short, 4>& polynomial) {
        for (polynomial_t<short, 4>::const_iterator it = polynomial.begin();
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

////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const exponent_t<float, 4>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<float, 4>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<float, 4>& polynomial) {
        for (polynomial_t<float, 4>::const_iterator it = polynomial.begin();
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

////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const exponent_t<double, 4>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<double, 4>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<double, 4>& polynomial) {
        for (polynomial_t<double, 4>::const_iterator it = polynomial.begin();
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

////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const polynomial_term_t<short, 4, mutation_tree_t>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<short, 4, mutation_tree_t>& polynomial) {
        for (polynomial_t<short, 4, mutation_tree_t>::const_iterator it = polynomial.begin();
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

size_t hash_value(const exponent_t<short, 2>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;

        return seed;
}

size_t hash_value(const exponent_t<short, 3>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;

        return seed;
}

size_t hash_value(const exponent_t<short, 4>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;

        return seed;
}


size_t hash_value(const exponent_t<short, 5>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;

        return seed;
}

size_t hash_value(const exponent_t<short, 10>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;
        seed += (size_t)exponent[5] << 10;
        seed += (size_t)exponent[6] << 12;
        seed += (size_t)exponent[7] << 14;
        seed += (size_t)exponent[8] << 16;
        seed += (size_t)exponent[9] << 18;

        return seed;
}

size_t hash_value(const exponent_t<float, 2>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;

        return seed;
}

size_t hash_value(const exponent_t<float, 3>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;

        return seed;
}

size_t hash_value(const exponent_t<float, 4>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;

        return seed;
}

size_t hash_value(const exponent_t<float, 5>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;

        return seed;
}

size_t hash_value(const exponent_t<float, 10>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;
        seed += (size_t)exponent[5] << 10;
        seed += (size_t)exponent[6] << 12;
        seed += (size_t)exponent[7] << 14;
        seed += (size_t)exponent[8] << 16;
        seed += (size_t)exponent[9] << 18;

        return seed;
}

size_t hash_value(const exponent_t<double, 2>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;

        return seed;
}

size_t hash_value(const exponent_t<double, 3>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;

        return seed;
}

size_t hash_value(const exponent_t<double, 4>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;

        return seed;
}

size_t hash_value(const exponent_t<double, 5>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;

        return seed;
}

size_t hash_value(const exponent_t<double, 10>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;
        seed += (size_t)exponent[5] << 10;
        seed += (size_t)exponent[6] << 12;
        seed += (size_t)exponent[7] << 14;
        seed += (size_t)exponent[8] << 16;
        seed += (size_t)exponent[9] << 18;

        return seed;
}
