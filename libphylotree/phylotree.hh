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

#ifndef PHYLOTREE_H
#define PHYLOTREE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cstring>
#include <boost/array.hpp>

#include <tfbayes/polynomial.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_node_t {
public:
        pt_node_t(short x = -1, double d = 0.0,
                  pt_node_t<CODE_TYPE, ALPHABET_SIZE>* left  = NULL,
                  pt_node_t<CODE_TYPE, ALPHABET_SIZE>* right = NULL) {
                this->x = x;
                this->d = d;
                this->left  = left;
                this->right = right;

                for (size_t i = 0; i < ALPHABET_SIZE; i++)
                        applicable[i] = false;
                applicable[ALPHABET_SIZE] = true;
        }

        bool leaf() const { return left == NULL && right == NULL; }
        bool root() const { return d == 0.0; }

        /* coded nucleotide */
        short x;
        /* distance to ancestor */
        double d;
        /* left child */
        pt_node_t<CODE_TYPE, ALPHABET_SIZE>* left;
        /* right child */
        pt_node_t<CODE_TYPE, ALPHABET_SIZE>* right;
        /* array that tells if it is possible to carry
         * a certain nucleotide upwards the tree */
        boost::array<bool, ALPHABET_SIZE+1> applicable;
        /* polynomial */
        boost::array<polynomial_t<CODE_TYPE, ALPHABET_SIZE>, ALPHABET_SIZE+1> carry;
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum;
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_root_t : public pt_node_t<CODE_TYPE, ALPHABET_SIZE> {
public:
        pt_root_t(short x = -1,
                  pt_node_t<CODE_TYPE, ALPHABET_SIZE>* left = NULL,
                  pt_node_t<CODE_TYPE, ALPHABET_SIZE>* right = NULL)
                : pt_node_t<CODE_TYPE, ALPHABET_SIZE>(x, 0.0, left, right) { }
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_leaf_t : public pt_node_t<CODE_TYPE, ALPHABET_SIZE> {
public:
        pt_leaf_t(short x, double d)
                : pt_node_t<CODE_TYPE, ALPHABET_SIZE>(x, d, NULL, NULL) { }
};

#define alphabet_size 4
typedef unsigned short code_t;

#endif /* PHYLOTREE_H */
