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

#include <Python.h>

#include <locale>
#include <cctype>
#include <sstream>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/format.hpp>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/operators.hpp>

#include <tfbayes/fg/distribution.hh>
#include <tfbayes/fg/variational.hh>
#include <tfbayes/fg/factor-graph.hh>
#include <tfbayes/interface/exceptions.hh>

using namespace boost::python;

// utility
// -----------------------------------------------------------------------------

static
std::vector<double> distribution_moment(const distribution_i& distribution, size_t n)
{
        switch (n) {
        case 1:
                return distribution.moment<1>();
        case 2:
                return distribution.moment<2>();
        default:
                return std::vector<double>();
        }
}

static
exponential_family_i* cast_exponential_family(const distribution_i& distribution)
{
        const exponential_family_i* result;

        if ((result = dynamic_cast<const exponential_family_i*>(&distribution))) {
                return result->clone();
        }
        raise_IOError("Distribution could not be casted to an exponential family!");
        // never reached
        return NULL;
}

static
factor_graph_t* construct_factor_graph(list& factor_nodes, list& variable_nodes, size_t threads = 1)
{
        factor_set_t fnodes;
        variable_set_t vnodes;

        for (ssize_t i = 0; i < len(factor_nodes); i++) {
                factor_node_i& ref = extract<factor_node_i&>(factor_nodes[i]);
                fnodes += &ref;
        }
        for (ssize_t i = 0; i < len(variable_nodes); i++) {
                variable_node_i& ref = extract<variable_node_i&>(variable_nodes[i]);
                vnodes += &ref;
        }
        factor_graph_t* fg = new factor_graph_t(fnodes, vnodes, threads);

        // make sure that no nodes are deleted
        fnodes.release().release();
        vnodes.release().release();

        return fg;
}
static
factor_graph_t* construct_factor_graph_1(list& factor_nodes, list& variable_nodes)
{
        return construct_factor_graph(factor_nodes, variable_nodes);
}
static
factor_graph_t* construct_factor_graph_2(list& factor_nodes, list& variable_nodes, size_t threads)
{
        return construct_factor_graph(factor_nodes, variable_nodes, threads);
}

static
const distribution_i& access_distribution(const factor_graph_t& factor_graph, const std::string& name)
{
        boost::optional<const distribution_i&> result = factor_graph[name];

        if (!result) {
                raise_IOError(boost::str(boost::format("Variable node `%s' not found in factor graph.")
                                         % name));
        }
        return *result;
}

static
void call_factor_graph(factor_graph_t& factor_graph)
{
        factor_graph();
}

static
void call_factor_graph(factor_graph_t& factor_graph, size_t n)
{
        factor_graph(n);
}

// interface
// -----------------------------------------------------------------------------

BOOST_PYTHON_MODULE(interface)
{
        // distributions
        // ---------------------------------------------------------------------
        class_<distribution_i, boost::noncopyable>("distribution_i", no_init)
                .def("moment", &distribution_moment)
                ;
        class_<dirac_distribution_t, bases<distribution_i> >("dirac_distribution_t", no_init)
                .def(init<std::vector<double> >())
                ;
        class_<exponential_family_i, bases<distribution_i>, boost::noncopyable>("exponential_family_i", no_init)
                .def("__init__",      make_constructor(&cast_exponential_family))
                .def("density",       &exponential_family_i::density)
                .def("base_measure",  &exponential_family_i::base_measure)
                .def("log_partition", &exponential_family_i::log_partition)
                .def("renormalize",   &exponential_family_i::renormalize)
                .def(self *= self)
                ;
        class_<normal_distribution_t, bases<exponential_family_i> >("normal_distribution_t")
                .def(init<double, double>())
                ;
        class_<gamma_distribution_t, bases<exponential_family_i> >("gamma_distribution_t")
                .def(init<double, double>())
                ;
        // factor nodes
        // ---------------------------------------------------------------------
        class_<factor_node_i, boost::noncopyable>("factor_node_i", no_init)
                .def("link", static_cast<bool (factor_node_i::*)(const std::string&, variable_node_i&)>(&factor_node_i::link))
                ;
        class_<normal_fnode_t, bases<factor_node_i> >("normal_fnode_t", no_init)
                .def(init<double, double>())
                .def(init<double, double, std::string>())
                ;
        class_<gamma_fnode_t, bases<factor_node_i> >("gamma_fnode_t", no_init)
                .def(init<double, double>())
                .def(init<double, double, std::string>())
                ;
        // variable nodes
        // ---------------------------------------------------------------------
        class_<variable_node_i, boost::noncopyable>("variable_node_i", no_init)
                .def("__call__", &variable_node_i::operator(), return_internal_reference<>())
                ;
        class_<data_vnode_t, bases<variable_node_i> >("data_vnode_t", no_init)
                .def(init<std::vector<double> >())
                .def(init<std::vector<double>, std::string>())
                ;
        class_<normal_vnode_t, bases<variable_node_i> >("normal_vnode_t")
                .def(init<std::string>())
                ;
        class_<gamma_vnode_t, bases<variable_node_i> >("gamma_vnode_t")
                .def(init<std::string>())
                ;
        // factor graph
        // ---------------------------------------------------------------------
        class_<factor_graph_t>("factor_graph_t", no_init)
                .def("__init__", make_constructor(&construct_factor_graph_1))
                .def("__init__", make_constructor(&construct_factor_graph_2))
                .def("__call__", static_cast<void (*)(factor_graph_t&)>(&call_factor_graph))
                .def("__call__", static_cast<void (*)(factor_graph_t&, size_t)>(&call_factor_graph))
                .def("__getitem__", &access_distribution, return_internal_reference<>())
                ;
}
