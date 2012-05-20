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

#include <phylotree-gradient-coefficient.hh>

using namespace std;

mutation_coefficient_t operator*(const mutation_coefficient_t& c1, const mutation_coefficient_t& c2) {
        mutation_coefficient_t tmp(c1);
        tmp *= c2;
        return tmp;
}

size_t
hash_value(const pmut_t& pmut)
{
        return (size_t)pmut.node;
}

ostream& operator<< (ostream& o, const mutation_product_t& product) {
        for (mutation_product_t::const_iterator it = product.begin(); it != product.end(); it++) {
                const pmut_t& mutation = *it;
                if (mutation) {
                        o << "M(" << mutation.node->name << ") ";
                }
                else {
                        o << "(1-M(" << mutation.node->name << ")) ";
                }
        }
        return o;
}

ostream& operator<< (ostream& o, const mutation_coefficient_t& coefficient) {

        o << "[ ";
        for (mutation_coefficient_t::const_iterator it = coefficient.begin(); it != coefficient.end(); it++) {
                if (it != coefficient.begin()) {
                        o << "+ ";
                }
                if (it->second != 1.0) {
                        o << it->second << " ";
                }
                o << it->first;
        }
        o << "]";
        return o;
}
