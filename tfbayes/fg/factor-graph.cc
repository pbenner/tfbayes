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
#include <iterator>

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
                while (true) {
                        // stop if maximum number of iterations is reached
                        if (n && (*n) == 0) break;
                        // decrement n
                        if (n) (*n)--;
                        // lock queue
                        mtx().lock();
                        // receive objects if queue is not empty
                        if (queue().empty()) {
                                mtx().unlock();
                                break;
                        }
                        else {
                                job = queue().front();
                                queue().pop();
                                mtx().unlock();
                                job->send_messages();
                        }
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

factor_graph_t::factor_graph_t(size_t threads) :
        _factor_nodes   (),
        _variable_nodes (),
        _threads        (threads)
{
}

factor_graph_t::factor_graph_t(const factor_set_t& factor_nodes,
                               const variable_set_t& variable_nodes,
                               size_t threads) :
        _factor_nodes   (),
        _variable_nodes (),
        _threads        (threads)
{
        // copy connectivity before proceeding
        clone_nodes(factor_nodes, variable_nodes);
}

factor_graph_t::factor_graph_t(const factor_graph_t& factor_graph) :
        _factor_nodes   (),
        _variable_nodes (),
        _queue          (),
        _threads        (factor_graph._threads)
{
        // clone all nodes of the factor graph
        clone_nodes(factor_graph._factor_nodes,
                    factor_graph._variable_nodes);
}

factor_graph_t&
factor_graph_t::operator+=(const factor_node_i& factor_node)
{
        return operator+=(factor_node.clone());
}

factor_graph_t&
factor_graph_t::operator+=(factor_node_i* factor_node)
{
        void (factor_graph_t::*tmp) (factor_graph_node_i*) = &factor_graph_t::add_node;

        _factor_nodes += factor_node;
        factor_node->observe(boost::bind(tmp, this, factor_node));

        return *this;
}

factor_graph_t&
factor_graph_t::operator+=(const variable_node_i& variable_node)
{
        return operator+=(variable_node.clone());
}

factor_graph_t&
factor_graph_t::operator+=(variable_node_i* variable_node)
{
        void (factor_graph_t::*tmp) (factor_graph_node_i*) = &factor_graph_t::add_node;

        _variable_nodes += variable_node;
        variable_node->observe(boost::bind(tmp, this, variable_node));

        return *this;
}

factor_graph_t&
factor_graph_t::operator+=(const factor_graph_t& factor_graph)
{
        clone_nodes(factor_graph._factor_nodes,
                    factor_graph._variable_nodes);

        return *this;
}

factor_graph_t&
factor_graph_t::replicate(size_t n)
{
        factor_graph_t tmp(*this);

        for (size_t i = 0; i < n; i++) {
                operator+=(tmp);
        }
        return *this;
}

bool
factor_graph_t::link(const std::string& fname, const std::string& which, const std::string& vname)
{
        bool result = false;

        for (factor_set_t::iterator it = _factor_nodes[fname];
             it != _factor_nodes.end() && it->name() == fname; it++) {
                for (variable_set_t::iterator is = _variable_nodes[vname];
                     is != _variable_nodes.end() && is->name() == vname; is++) {
                        result |= it->link(which, *is);
                }
        }
        return result;
}

void
factor_graph_t::operator()(boost::optional<size_t> n) {
        vector<boost::thread*> threads(_threads);

        // initialize the network by letting all variable nodes send
        // their messages first
        for (variable_set_t::iterator it = _variable_nodes.begin();
             it != _variable_nodes.end(); it++) {
                it->send_messages();
        }

        // sample
        for (size_t i = 0; i < _threads; i++) {
                threads[i] = new boost::thread(factor_graph_thread_t(_queue, mtx), n);
        }
        // join threads
        for (size_t i = 0; i < _threads; i++) {
                threads[i]->join();
        }
        for (size_t i = 0; i < _threads; i++) {
                delete(threads[i]);
        }
}

boost::optional<const distribution_i&>
factor_graph_t::distribution(const string& name, size_t i) const
{
        variable_set_t::const_iterator it = _variable_nodes[name];
        // advance iterator i steps
        advance(it, i);

        if (it == _variable_nodes.cend() ||
            it->name() != name) {
                return boost::optional<const distribution_i&>();
        }
        return (*it)();
}

boost::optional<const variable_node_i&>
factor_graph_t::variable_node(const string& name, size_t i) const
{
        variable_set_t::const_iterator it = _variable_nodes[name];
        // advance iterator i steps
        advance(it, i);

        if (it == _variable_nodes.cend() ||
            it->name() != name) {
                return boost::optional<const variable_node_i&>();
        }
        return *it;
}

boost::optional<variable_node_i&>
factor_graph_t::variable_node(const string& name, size_t i)
{
        variable_set_t::iterator it = _variable_nodes[name];
        // advance iterator i steps
        advance(it, i);

        if (it == _variable_nodes.end() ||
            it->name() != name) {
                return boost::optional<variable_node_i&>();
        }
        return *it;
}

boost::optional<const data_vnode_t&>
factor_graph_t::data_vnode(const string& name, size_t i) const
{
        variable_set_t::const_iterator it = _variable_nodes[name];
        // advance iterator i steps
        advance(it, i);

        if (it == _variable_nodes.cend() ||
            it->name() != name) {
                return boost::optional<const data_vnode_t&>();
        }
        return static_cast<const data_vnode_t&>(*it);
}

boost::optional<data_vnode_t&>
factor_graph_t::data_vnode(const string& name, size_t i)
{
        variable_set_t::iterator it = _variable_nodes[name];
        // advance iterator i steps
        advance(it, i);

        if (it == _variable_nodes.end() ||
            it->name() != name          ||
            typeid(*it) != typeid(data_vnode_t)) {
                return boost::optional<data_vnode_t&>();
        }
        return static_cast<data_vnode_t&>(*it);
}

void
factor_graph_t::add_node(factor_graph_node_i* node) {
        mtx.lock();
        debug("adding node " << typeid(*node).name() << endl);
        _queue.push(node);
        mtx.unlock();
}

void
factor_graph_t::clone_nodes(const factor_set_t& fnodes,
                            const variable_set_t& vnodes)
{
        // map old nodes to new ones
        std::map<const variable_node_i*, variable_node_i*> vmap;
        std::map<const factor_node_i*,   factor_node_i*>   fmap;

        // clone variable nodes
        for (variable_set_t::const_iterator it = vnodes.cbegin();
             it != vnodes.cend(); it++) {
                variable_node_i* node = it->clone();
                operator+=(node);
                // keep track of which node replacements
                vmap[&*it] = node;
        }
        // clone factor nodes
        for (factor_set_t::const_iterator it = fnodes.cbegin();
             it != fnodes.cend(); it++) {
                factor_node_i* node = it->clone();
                operator+=(node);
                // keep track of which node replacements
                fmap[&*it] = node;
        }
        // copy connectivity
        for (factor_set_t::const_iterator it = fnodes.cbegin();
             it != fnodes.cend(); it++) {
                const vector<variable_node_i*>& neighbors =
                        it->neighbors();
                for (size_t j = 0; j < neighbors.size(); j++) {
                        if (!neighbors[j])
                                continue;
                        fmap[&*it]->link(j, *vmap[neighbors[j]]);
                }
        }
}
