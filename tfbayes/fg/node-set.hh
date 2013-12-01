/* Copyright (C) 2013 Philipp Benner
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

#ifndef __TFBAYES_FG_NODE_SET_HH__
#define __TFBAYES_FG_NODE_SET_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <string>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <tfbayes/fg/node-types.hh>
#include <tfbayes/utility/default-operator.hh>

template<typename T1, typename T2>
T2 get_value(T1 pair) {
        return *pair.second;
}


template <class T>
class node_set_t : public boost::ptr_multimap<std::string, T> {
public:
        // type of the base class
        typedef boost::ptr_multimap<std::string, T> base_t;

        // type of the transformed iterator
        typedef boost::transform_iterator<
                const T& (*)(typename base_t::const_reference),
                typename base_t::const_iterator
                > const_iterator;
        typedef boost::transform_iterator<
                T& (*)(typename base_t::reference),
                typename base_t::iterator
                > iterator;

        friend void swap(node_set_t<T>& left, node_set_t<T>& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
        }
        default_assignment_operator(node_set_t)

        // access values
        const_iterator operator[](const std::string& name) const {
                return boost::make_transform_iterator(base_t::find(name), get_value<typename base_t::const_reference, const T&>);
        }
        iterator operator[](const std::string& name) {
                return boost::make_transform_iterator(base_t::find(name), get_value<typename base_t::reference, T&>);
        }
        const_iterator cbegin() const {
                return boost::make_transform_iterator(base_t::begin(), get_value<typename base_t::const_reference, const T&>);
        }
        iterator begin() {
                return boost::make_transform_iterator(base_t::begin(), get_value<typename base_t::reference, T&>);
        }
        const_iterator cend() const {
                return boost::make_transform_iterator(base_t::end(), get_value<typename base_t::const_reference, const T&>);
        }
        iterator end() {
                return boost::make_transform_iterator(base_t::end(), get_value<typename base_t::reference, T&>);
        }
        node_set_t<T>& operator+=(T* node) {
                std::string tmp(node->name());
                base_t::insert(tmp, node);
                return *this;
        }
};

typedef node_set_t<factor_node_i> factor_set_t;
typedef node_set_t<variable_node_i> variable_set_t;

#endif /* __TFBAYES_FG_NODE_SET_HH__ */
