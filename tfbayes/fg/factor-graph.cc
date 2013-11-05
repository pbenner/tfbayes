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

#include <iostream>

#include <boost/bind.hpp>

#include <factor-graph.hh>

// class that implements threads on factor graphs
////////////////////////////////////////////////////////////////////////////////

class factor_graph_thread_t {
public:
        factor_graph_thread_t(std::queue<factor_graph_node_i*>& queue,
                              boost::mutex& mtx) :
                _queue(&queue),
                _mtx  (&mtx)
                { }

        void operator()(boost::optional<size_t> n = boost::optional<size_t>()) {
                // pointer to a node
                factor_graph_node_i* job;
                // get jobs from the queue
                while ( !queue().empty()) {
                        // stop if maximum number of iterations is reached
                        if (n && (*n) == 0) break;
                        // decrement n
                        if (n) (*n)--;
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
        boost::mutex& mtx() {
                return *_mtx;
        }

protected:
        std::queue<factor_graph_node_i*>* _queue;
        boost::mutex* _mtx;
};

// the factor graph
////////////////////////////////////////////////////////////////////////////////

using namespace std;

factor_graph_t::factor_graph_t(const vector<variable_node_i*> variable_nodes,
                               const vector<  factor_node_i*> factor_nodes,
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
        _variable_nodes (),
        _factor_nodes   (),
        _queue          (),
        _threads        (factor_graph._threads)
{
        std::map<variable_node_i*, variable_node_i*> vmap;

        // clone all nodes of the factor graph
        for (size_t i = 0; i < factor_graph._variable_nodes.size(); i++) {
                _variable_nodes.push_back(factor_graph._variable_nodes[i]->clone());
                // keep track of which node replacements
                vmap[factor_graph._variable_nodes[i]] = _variable_nodes[i];
        }
        for (size_t i = 0; i < factor_graph._factor_nodes.size(); i++) {
                _factor_nodes.push_back(factor_graph._factor_nodes[i]->clone());
        }
        // find connected nodes and link them
        for (size_t i = 0; i < factor_graph._factor_nodes.size(); i++) {
                const vector<variable_node_i*>& neighbors =
                        factor_graph._factor_nodes[i]->neighbors();
                for (size_t j = 0; j < neighbors.size(); j++) {
                        if (!neighbors[j])
                                continue;
                        _factor_nodes[i]->link(j, *vmap[neighbors[j]]);
                }
        }
}

factor_graph_t&
factor_graph_t::operator=(const factor_graph_t& node) {
        factor_graph_t tmp(node);
        swap(*this, tmp);
        return *this;
}

void
factor_graph_t::operator()(boost::optional<size_t> n) {
        vector<boost::thread*> threads(_threads);

        cout << "creating " << _threads << " threads" << endl;
        cout << "n iterations " << n << endl;
        cout << "n iterations " << std::numeric_limits<size_t>::infinity() << endl;

        // add all nodes to the queue
        for (size_t i = 0; i < _variable_nodes.size(); i++) {
                _queue.push(_variable_nodes[i]);
        }
        for (size_t i = 0; i < _factor_nodes.size(); i++) {
                _queue.push(_factor_nodes[i]);
        }

        // sample
        for (size_t i = 0; i < _threads; i++) {
                if (n) {
                        threads[i] = new boost::thread(factor_graph_thread_t(_queue, mtx), (*n)/_threads);
                }
                else {
                        threads[i] = new boost::thread(factor_graph_thread_t(_queue, mtx));
                }
        }
        // join threads
        for (size_t i = 0; i < _threads; i++) {
                threads[i]->join();
        }
        for (size_t i = 0; i < _threads; i++) {
                delete(threads[i]);
        }
}
