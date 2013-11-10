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

#include <map>

#include <boost/optional.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include <tfbayes/fg/node-types.hh>

template <typename T>
class node_set_t : public boost::ptr_vector<T> {
public:
        typedef boost::ptr_vector<T> base_t;
        typedef std::map<std::string, size_t> map_t;

        node_set_t<T>& operator+=(T* node) {
                base_t::push_back(node);
                // insert node into the map if a name is available
                if (node->name() != "") {
                        _map[node->name()] = base_t::size()-1;
                }
                return *this;
        }
        using base_t::operator[];
        boost::optional<T&> operator[](const std::string& name) {
                map_t::iterator it = _map.find(name);
                if (it != _map.end()) {
                        return base_t::operator[](it->second);
                }
                return boost::optional<T&>();
        }
        boost::optional<const T&> operator[](const std::string& name) const {
                map_t::const_iterator it = _map.find(name);
                if (it != _map.end()) {
                        return base_t::operator[](it->second);
                }
                return boost::optional<const T&>();
        }

protected:
        map_t _map;
};

typedef node_set_t<factor_node_i> factor_set_t;
typedef node_set_t<variable_node_i> variable_set_t;

#endif /* __TFBAYES_FG_NODE_SET_HH__ */
