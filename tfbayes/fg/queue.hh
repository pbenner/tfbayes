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

#ifndef __TFBAYES_FG_QUEUE_HH__
#define __TFBAYES_FG_QUEUE_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <set>

#include <boost/thread.hpp>

#include <tfbayes/fg/node-types.hh>

class fg_queue_t : public boost::mutex, public std::set<factor_graph_node_i*> {
public:
        typedef std::set<factor_graph_node_i*> base_t;

        void push(factor_graph_node_i* node) {
                base_t::insert(node);
        }
        factor_graph_node_i* pop() {
                factor_graph_node_i* node = *base_t::begin();
                base_t::erase(node);
                return node;
        }
};

#endif /* __TFBAYES_FG_QUEUE_HH__ */
