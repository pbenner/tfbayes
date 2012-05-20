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

#ifndef PHYLOTREE_GRADIENT_COEFFICIENT_HH
#define PHYLOTREE_GRADIENT_COEFFICIENT_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <set>
#include <boost/unordered_map.hpp>
#include <ostream>

#include <phylotree.hh>

class pmut_t {
public:
        pmut_t(pt_node_t* node, bool mutation = true)
                : node(node), mutation(mutation) { }

        operator bool() const {
                return mutation;
        }
        bool operator==(const pmut_t& m) const {
                return node == m.node && mutation == m.mutation;
        }
        bool operator<(const pmut_t& m) const {
                if (node < m.node) {
                        return true;
                }
                if (node == m.node) {
                        return mutation < m.mutation;
                }
                return false;
        }

        pt_node_t* node;

private:
        bool mutation;
};

class mutation_product_t : public std::multiset<pmut_t> {
public:
        mutation_product_t()
                : std::multiset<pmut_t>() {}
        mutation_product_t(const pmut_t& mutation)
                : std::multiset<pmut_t>() {
                insert(mutation);
        }

        mutation_product_t operator*=(const pmut_t& mutation) {
                insert(mutation);
                return *this;
        }
        mutation_product_t operator*=(const mutation_product_t& product) {

                for (mutation_product_t::const_iterator it = product.begin(); it != product.end(); it++) {
                        operator*=(*it);
                }
                return *this;
        }
};

class mutation_coefficient_t : public boost::unordered_map<mutation_product_t, double> {
public:
        mutation_coefficient_t()
                : boost::unordered_map<mutation_product_t, double>() {}
        mutation_coefficient_t(double constant)
                : boost::unordered_map<mutation_product_t, double>() {
                operator[](mutation_product_t()) += constant;
        }
        mutation_coefficient_t(const pmut_t& mutation)
                : boost::unordered_map<mutation_product_t, double>() {
                operator[](mutation) += 1.0;
        }
        mutation_coefficient_t(const mutation_product_t& product)
                : boost::unordered_map<mutation_product_t, double>() {
                operator[](product) += 1.0;
        }

        operator bool() const {
                return size();
        }
        mutation_coefficient_t operator-() const {
                mutation_coefficient_t coefficient;
                for (mutation_coefficient_t::const_iterator it = begin(); it != end(); it++) {
                        coefficient[it->first] = -(find(it->first)->second);
                }
                return coefficient;
        }
        mutation_coefficient_t operator+=(
                const mutation_coefficient_t& coefficient) {

                for (mutation_coefficient_t::const_iterator it = coefficient.begin(); it != coefficient.end(); it++) {
                        operator[](it->first) += it->second;
                        if (operator[](it->first) == 0.0) {
                                erase(it->first);
                        }
                }
                return *this;
        }
        mutation_coefficient_t operator*=(
                const mutation_coefficient_t& coefficient) {

                mutation_coefficient_t tmp;
                for (mutation_coefficient_t::const_iterator it = begin(); it != end(); it++) {
                        for (mutation_coefficient_t::const_iterator is = coefficient.begin(); is != coefficient.end(); is++) {
                                mutation_product_t product(it->first);
                                product *= (is->first);
                                tmp[product] += (it->second)*(is->second);
                        }
                }
                operator=(tmp);

                return *this;
        }
};

//mutation_coefficient_t operator*(const mutation_coefficient_t& c1, const mutation_coefficient_t& c2);

size_t hash_value(const pmut_t& pmut);

std::ostream& operator<< (std::ostream& o, const mutation_product_t& product);
std::ostream& operator<< (std::ostream& o, const mutation_coefficient_t& coefficient);

#endif /* PHYLOTREE_GRADIENT_COEFFICIENT_HH */
