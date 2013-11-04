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

#include <boost/bind.hpp>
#include <boost/thread.hpp>

#include <factor-graph.hh>

// class that implements threads on factor graphs
////////////////////////////////////////////////////////////////////////////////

class factor_graph_thread_t {
public:
        factor_graph_thread_t(std::queue<factor_graph_node_i*>& queue,
                              std::mutex& mtx) :
                _queue(&queue),
                _mtx  (&mtx)
                { }

        void operator()(size_t n) {
                // pointer to a node
                factor_graph_node_i* job;
                // get jobs from the queue
                for (; !queue().empty() && n > 0; n--) {
                        mtx().lock();
                        job = queue().front();
                        queue().pop();
                        mtx().unlock();
                        job->send_messages();
                }
        }
        std::queue<factor_graph_node_i*>& queue() {
                return *_queue;
        }
        std::mutex& mtx() {
                return *_mtx;
        }

protected:
        std::queue<factor_graph_node_i*>* _queue;
        std::mutex* _mtx;
};

// the factor graph
////////////////////////////////////////////////////////////////////////////////

factor_graph_t::factor_graph_t(const std::vector<variable_node_i*> variable_nodes,
               const std::vector<  factor_node_i*> factor_nodes,
               size_t threads) :
        _variable_nodes (variable_nodes),
        _factor_nodes   (factor_nodes),
        _threads        (threads)
{
        void (factor_graph_t::*tmp) (factor_graph_node_i*) = &factor_graph_t::add_node;

        for (size_t i = 0; i < _variable_nodes.size(); i++) {
                _variable_nodes[i]->observe(
                        boost::bind(tmp, this, _variable_nodes[i]));
        }
        for (size_t i = 0; i < _factor_nodes.size(); i++) {
                _factor_nodes[i]->observe(
                        boost::bind(tmp, this, _factor_nodes[i]));
        }
}

factor_graph_t::~factor_graph_t() {
        for (size_t i = 0; i < _variable_nodes.size(); i++) {
                delete(_variable_nodes[i]);
        }
        for (size_t i = 0; i < _factor_nodes.size(); i++) {
                delete(_factor_nodes[i]);
        }
}

factor_graph_t::factor_graph_t(const factor_graph_t& factor_graph) :
        _variable_nodes (factor_graph._variable_nodes),
        _factor_nodes   (factor_graph._factor_nodes),
        _queue          (factor_graph._queue),
        _threads        (factor_graph._threads)
{ }

factor_graph_t&
factor_graph_t::operator=(const factor_graph_t& node) {
        using std::swap;
        factor_graph_t tmp(node);
        swap(*this, tmp);
        return *this;
}

void
factor_graph_t::operator()(size_t n) {
        std::vector<boost::thread*> threads(_threads);

        // add all nodes to the queue
        for (size_t i = 0; i < _variable_nodes.size(); i++) {
                _queue.push(_variable_nodes[i]);
        }
        for (size_t i = 0; i < _factor_nodes.size(); i++) {
                _queue.push(_factor_nodes[i]);
        }

        // sample
        for (size_t i = 0; i < _threads; i++) {
                threads[i] = new boost::thread(factor_graph_thread_t(_queue, mtx), n/_threads);
        }
        // join threads
        for (size_t i = 0; i < _threads; i++) {
                threads[i]->join();
        }
        for (size_t i = 0; i < _threads; i++) {
                delete(threads[i]);
        }
}
