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

#include <boost/ptr_container/ptr_vector.hpp>

#include <tfbayes/fg/node-types.hh>

template <typename T>
class node_set_t : public boost::ptr_vector<T> {
public:
        node_set_t<T>& operator+=(T* node) {
                this->push_back(node);
                return *this;
        }
};

typedef node_set_t<factor_node_i> factor_set_t;
typedef node_set_t<variable_node_i> variable_set_t;

#endif /* __TFBAYES_FG_NODE_SET_HH__ */
