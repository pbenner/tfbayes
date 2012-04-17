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

#include <incomplete-expression.hh>

using namespace std;

ostream& operator<< (ostream& o, const incomplete_nodeset_t& exponent) {
        for (incomplete_nodeset_t::const_iterator it = exponent.begin(); it != exponent.end(); it++) {
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

ostream& operator<< (ostream& o, const incomplete_nodeterm_t& term) {
        if (term.coefficient() != 1.0) {
                o << term.coefficient()
                  << " ";
        }
        o << (const incomplete_nodeset_t&)term;

        return o;
}

ostream& operator<< (ostream& o, const incomplete_expression_t& expression) {
        for (incomplete_expression_t::const_iterator it = expression.begin(); it != expression.end(); it++) {
                if (it != expression.begin()) {
                        o << " + ";
                }
                o << *it;
        }

        return o;
}

size_t hash_value(const incomplete_nodeset_t& nodeset) {
        size_t seed = 0;
        for (incomplete_nodeset_t::const_iterator it = nodeset.begin(); it != nodeset.end(); it++) {
                boost::hash_combine(seed, hash_value(*it));
        }
        boost::hash_combine(seed, hash_value(nodeset.incomplete()));
        return seed;
}

incomplete_nodeterm_t operator*(double constant, const incomplete_nodeterm_t& term) {
        incomplete_nodeterm_t result(term);
        result *= constant;
        return result;
}
incomplete_nodeterm_t operator*(const incomplete_nodeterm_t& term, double constant) {
        incomplete_nodeterm_t result(term);
        result *= constant;
        return result;
}
incomplete_nodeterm_t operator*(const incomplete_nodeterm_t& term1, const incomplete_nodeterm_t& term2) {
        incomplete_nodeterm_t result(term1);
        result *= term2;
        return result;
}
incomplete_expression_t operator*(const incomplete_expression_t& expression1, const incomplete_expression_t& expression2) {
        incomplete_expression_t result(expression1);
        result *= expression2;
        return result;
}
