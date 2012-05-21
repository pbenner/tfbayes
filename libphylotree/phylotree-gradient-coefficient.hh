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

#include <phylotree.hh>

class pmut_t {
public:
        pmut_t()
                : node(NULL), mutation(false) { }
        pmut_t(pt_node_t* node, bool mutation = true)
                : node(node), mutation(mutation) { }

        bool operator!() const {
                return node == NULL;
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
        double eval() const {
                if (operator!()) {
                        return 1.0;
                }
                if (mutation) {
                        return 1.0-exp(-node->d);
                }
                else {
                        return exp(-node->d);
                }
        }

        pt_node_t* node;
        bool mutation;
};

class mutation_tree_t {
public:
        mutation_tree_t()
                : _mult(false), _sum(false), _pmut(), _constant(0.0),
                  _left(NULL ), _right(NULL) { }
        mutation_tree_t(double c)
                : _mult(false), _sum(false), _pmut(), _constant(c),
                  _left(NULL ), _right(NULL) { }
        mutation_tree_t(const mutation_tree_t& tree)
                : _mult(tree._mult), _sum(tree._sum), _pmut(tree._pmut), _constant(tree._constant),
                  _left(NULL), _right(NULL) {
                if (tree._left ) { _left  = new mutation_tree_t(*tree._left ); }
                if (tree._right) { _right = new mutation_tree_t(*tree._right); }
        }
        mutation_tree_t(const pmut_t& pmut)
                : _mult(false), _sum(false), _pmut(pmut), _constant(1.0),
                  _left(NULL), _right(NULL) { }
        ~mutation_tree_t() {
                if (!_left ) { delete(_left ); _left  = NULL; }
                if (!_right) { delete(_right); _right = NULL; }
        }

        bool operator==(const mutation_tree_t& tree) const {
                bool result = true;
                if ((_left  && tree._left  == NULL) || (_left  == NULL && tree._left )) {
                        return false;
                }
                if ((_right && tree._right == NULL) || (_right == NULL && tree._right)) {
                        return false;
                }
                result = result && (this->sum     () == tree.sum     ());
                result = result && (this->mult    () == tree.mult    ());
                result = result && (this->pmut    () == tree.pmut    ());
                result = result && (this->constant() == tree.constant());
                if (_left  && tree._left ) {
                        result = result && (this->left()  == tree.left ());
                }
                if (_right && tree._right) {
                        result = result && (this->right() == tree.right());
                }
                return result;
        }
        bool operator!=(const mutation_tree_t& tree) const {
                return !operator==(tree);
        }
        bool operator!() const {
                return constant() == 0.0;
        }
        mutation_tree_t operator-() const {

                mutation_tree_t tree(*this);
                tree *= -1.0;
                return tree;
        }
        mutation_tree_t& operator*=(double c) {
                _constant *= c;
                return *this;
        }
        mutation_tree_t& operator*=(const mutation_tree_t& tree) {
                if (tree  == mutation_tree_t(1.0)) {
                        return *this;
                }
                if (*this == mutation_tree_t(1.0)) {
                        *this = tree;
                        return *this;
                }

                mutation_tree_t* copy_left  = new mutation_tree_t(*this);
                mutation_tree_t* copy_right = new mutation_tree_t(tree);

                _mult     = true;
                _sum      = false;
                _pmut     = pmut_t();
                _constant = 1.0;
                _left     = copy_left;
                _right    = copy_right;

                return *this;
        }
        mutation_tree_t& operator+=(const mutation_tree_t& tree) {
                if ( tree == mutation_tree_t(0.0)) {
                        return *this;
                }
                if (*this == mutation_tree_t(0.0)) {
                        *this = tree;
                        return *this;
                }

                mutation_tree_t* copy_left  = new mutation_tree_t(*this);
                mutation_tree_t* copy_right = new mutation_tree_t(tree);

                _mult     = false;
                _sum      = true;
                _pmut     = pmut_t();
                _constant = 1.0;
                _left     = copy_left;
                _right    = copy_right;

                return *this;
        }

        bool mult() const {
                return _mult;
        }
        bool sum() const {
                return _sum;
        }
        bool leaf() const {
                return !(mult() || sum());
        }
        const pmut_t& pmut() const {
                return _pmut;
        }
        double& constant() {
                return _constant;
        }
        const double& constant() const {
                return _constant;
        }
        mutation_tree_t& left() {
                return *_left;
        }
        const mutation_tree_t& left() const {
                return *_left;
        }
        const mutation_tree_t& right() const {
                return *_right;
        }
        double eval() const {
                if (leaf()) {
                        if (!!pmut()) {
                                return constant()*pmut().eval();
                        }
                        else {
                                return pmut().eval();
                        }
                }
                else {
                        double result1 = left() .eval();
                        double result2 = right().eval();
                        if (mult()) {
                                return constant()*result1*result2;
                        }
                        else {
                                return constant()*(result1+result2);
                        }
                }
        }

private:
        bool _mult;
        bool _sum;
        pmut_t _pmut;
        double _constant;

        mutation_tree_t* _left;
        mutation_tree_t* _right;
};

size_t hash_value(const pmut_t& pmut);

#include <ostream>

std::ostream& operator<< (std::ostream& o, const mutation_tree_t& tree);

#endif /* PHYLOTREE_GRADIENT_COEFFICIENT_HH */
