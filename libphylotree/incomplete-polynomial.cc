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

#include <incomplete-polynomial.hh>

using namespace std;

ostream& operator<< (ostream& o, const node_set_t& node_set) {
        for (node_set_t::const_iterator it = node_set.begin(); it != node_set.end(); it++) {
                if (it != node_set.begin()) {
                        o << ",";
                }
                o << (*it)->name;
        }

        return o;
}

ostream& operator<< (ostream& o, const incomplete_exponent_t& exponent) {
        for (incomplete_exponent_t::const_iterator it = exponent.begin(); it != exponent.end(); it++) {
                o << "phi("
                  << *it
                  << ")";
        }
        if (!exponent.incomplete().empty()) {
                o << "I("
                  << exponent.incomplete()
                  << ")";
        }
        return o;
}

ostream& operator<< (ostream& o, const incomplete_term_t& term) {
        if (term.coefficient() != 1.0) {
                o << term.coefficient()
                  << " ";
        }
        o << (const incomplete_exponent_t&)term;

        return o;
}

ostream& operator<< (ostream& o, const incomplete_polynomial_t& polynomial) {
        for (incomplete_polynomial_t::const_iterator it = polynomial.begin(); it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + ";
                }
                o << *it;
        }

        return o;
}

size_t hash_value(const node_set_t& set) {
        size_t seed = 0;
        for (node_set_t::const_iterator it = set.begin(); it != set.end(); it++) {
                seed += (size_t)*it;
        }
        return seed;
}

size_t hash_value(const incomplete_exponent_t& exponent) {
        size_t seed = 0;
        for (incomplete_exponent_t::const_iterator it = exponent.begin(); it != exponent.end(); it++) {
                boost::hash_combine(seed, hash_value(*it));
        }
        boost::hash_combine(seed, hash_value(exponent.incomplete()));
        return seed;
}

incomplete_term_t operator*(double constant, const incomplete_term_t& term) {
        incomplete_term_t result(term);
        result *= constant;
        return result;
}
incomplete_term_t operator*(const incomplete_term_t& term, double constant) {
        incomplete_term_t result(term);
        result *= constant;
        return result;
}
incomplete_term_t operator*(const incomplete_term_t& term1, const incomplete_term_t& term2) {
        incomplete_term_t result(term1);
        result *= term2;
        return result;
}
incomplete_polynomial_t operator*(const incomplete_polynomial_t& poly1, const incomplete_polynomial_t& poly2) {
        incomplete_polynomial_t result(poly1);
        result *= poly2;
        return result;
}
