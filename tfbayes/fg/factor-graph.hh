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

#ifndef __TFBAYES_FG_FACTOR_GRAPH_HH__
#define __TFBAYES_FG_FACTOR_GRAPH_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <limits>
#include <queue>
#include <vector>

#include <boost/thread.hpp>

#include <node-types.hh>

class factor_graph_t {
public:
        factor_graph_t(const std::vector<variable_node_i*> variable_nodes,
                       const std::vector<  factor_node_i*> factor_nodes,
                       size_t threads = 1);
        factor_graph_t(const factor_graph_t& factor_graph);

        virtual ~factor_graph_t();


        virtual factor_graph_t* clone() const {
                return new factor_graph_t(*this);
        }

        friend void swap(factor_graph_t& left, factor_graph_t& right) {
                using std::swap;
                swap(left._variable_nodes, right._variable_nodes);
                swap(left._factor_nodes,   right._factor_nodes);
                swap(left._queue,          right._queue);
                swap(left._threads,        right._threads);
        }

        factor_graph_t& operator=(const factor_graph_t& node);

        // execute the message passing algorithm
        void operator()(size_t n = std::numeric_limits<size_t>::infinity());

protected:
        std::vector<variable_node_i*> _variable_nodes;
        std::vector<  factor_node_i*> _factor_nodes;
        // queue of nodes that need to send messages
        std::queue<factor_graph_node_i*> _queue;
        // number of threads
        size_t _threads;
private:
        // add a node to the queue
        void add_node(factor_graph_node_i* node) {
                mtx.lock();
                _queue.push(node);
                mtx.unlock();
        }
        // lock for the queue
        mutable boost::mutex mtx;
};

#endif /* __TFBAYES_FG_FACTOR_GRAPH_HH__ */
