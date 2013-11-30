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

#include <tfbayes/fg/factor-graph.hh>
#include <tfbayes/utility/strtools.hh>

// the factor graph
////////////////////////////////////////////////////////////////////////////////

using namespace std;

factor_graph_t::factor_graph_t() :
        _factor_nodes   (),
        _variable_nodes (),
        _factor_queue   (),
        _variable_queue (),
        _energy_cache   ()
{
}

factor_graph_t::factor_graph_t(const factor_graph_t& factor_graph) :
        _factor_nodes   (),
        _variable_nodes (),
        _factor_queue   (),
        _variable_queue (),
        _energy_cache   ()
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
        _factor_nodes += factor_node;
        factor_node->observe(boost::bind(&factor_graph_t::add_factor_node, this, factor_node));

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
        _variable_nodes += variable_node;
        variable_node->observe(boost::bind(&factor_graph_t::add_variable_node, this, variable_node));

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
factor_graph_t::link(const string& tag, const std::string& vname)
{
        if (token_first(tag, ':').size() < 2) {
                return false;
        }
        string tag1 = token_first(tag, ':')[0];
        string tag2 = token_first(tag, ':')[1];

        bool result = false;

        for (factor_set_t::iterator it = _factor_nodes[tag1];
             it != _factor_nodes.end() && it->name() == tag1; it++) {
                for (variable_set_t::iterator is = _variable_nodes[vname];
                     is != _variable_nodes.end() && is->name() == vname; is++) {
                        result |= it->link(tag2, *is);
                }
        }
        return result;
}

vector<double>
factor_graph_t::operator()(boost::optional<size_t> n) {
        // a cache for the free energy
        // and a thread safe vector for the results
        vector<double> history;
        // limit the number of jobs
        _variable_queue.set_limit(n);
        // initialize queue
        for (factor_set_t::iterator it = _factor_nodes.begin();
             it != _factor_nodes.end(); it++) {
                _energy_cache.update(&*it, it->free_energy());
        }
        for (variable_set_t::iterator it = _variable_nodes.begin();
             it != _variable_nodes.end(); it++) {
                it->notify();
                _energy_cache.update(&*it, it->free_energy());
        }
        history.push_back(free_energy());
        // iterate network
        while (variable_node_i* job = _variable_queue.pop()) {
                debug("--------------------------------------------------------------------------------"
                      << std::endl);
                job->update();
                history.push_back(free_energy());
        }
        return history;
}

double
factor_graph_t::free_energy()
{
        while (factor_node_i* job = _factor_queue.pop()) {
                _energy_cache.update(job, job->free_energy());
        }
        return _energy_cache();
}

boost::optional<const exponential_family_i&>
factor_graph_t::distribution(const string& name, size_t i) const
{
        variable_set_t::const_iterator it = _variable_nodes[name];
        // advance iterator i steps
        advance(it, i);

        if (it == _variable_nodes.cend() ||
            it->name() != name) {
                return boost::optional<const exponential_family_i&>();
        }
        return it->distribution();
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

void
factor_graph_t::add_factor_node(factor_node_i* factor_node) {
        debug(boost::format("*** adding factor node %s:%x to the queue ***\n")
              % factor_node->name() % factor_node);
        _factor_queue.push(factor_node);
}

void
factor_graph_t::add_variable_node(variable_node_i* variable_node) {
        debug(boost::format("*** adding variable node %s:%x to the queue ***\n")
              % variable_node->name() % variable_node);
        _variable_queue.push(variable_node);
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
                // keep track of node replacements
                vmap[&*it] = node;
        }
        // clone factor nodes
        for (factor_set_t::const_iterator it = fnodes.cbegin();
             it != fnodes.cend(); it++) {
                factor_node_i* node = it->clone();
                operator+=(node);
                // keep track of node replacements
                fmap[&*it] = node;
        }
        // copy connectivity
        for (factor_set_t::const_iterator it = fnodes.cbegin();
             it != fnodes.cend(); it++) {
                const factor_node_i::neighbors_t& neighbors =
                        it->neighbors();
                for (size_t j = 0; j < neighbors.size(); j++) {
                        if (!neighbors[j])
                                continue;
                        fmap[&*it]->link(neighbors[j], *vmap[neighbors[j]]);
                }
        }
}
