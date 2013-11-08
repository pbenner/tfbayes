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
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/operators.hpp>

#include <distribution.hh>
#include <variational.hh>
#include <factor-graph.hh>

using namespace boost::python;

// utility
// -----------------------------------------------------------------------------

static
double distribution_moment(const distribution_i& distribution, size_t n)
{
        switch (n) {
        case 1:
                return distribution.moment<1>();
        case 2:
                return distribution.moment<2>();
        default:
                return std::numeric_limits<double>::infinity();
        }
}

static
factor_graph_t* construct_factor_graph(list& factor_nodes, list& variable_nodes)
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
        return new factor_graph_t(fnodes, vnodes);
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
                .def(init<double>())
                ;
        class_<exponential_family_i, bases<distribution_i>, boost::noncopyable>("exponential_family_i", no_init)
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
                ;
        class_<gamma_fnode_t, bases<factor_node_i> >("gamma_fnode_t", no_init)
                .def(init<double, double>())
                ;
        // variable nodes
        // ---------------------------------------------------------------------
        class_<variable_node_i, boost::noncopyable>("variable_node_i", no_init)
                .def("__call__", &variable_node_i::operator(), return_internal_reference<>())
                ;
        class_<data_vnode_t, bases<variable_node_i> >("data_vnode_t", no_init)
                .def(init<double>())
                ;
        class_<normal_vnode_t, bases<variable_node_i> >("normal_vnode_t")
                ;
        class_<gamma_vnode_t, bases<variable_node_i> >("gamma_vnode_t")
                ;
        // factor graph
        // ---------------------------------------------------------------------
        class_<factor_graph_t>("factor_graph_t", no_init)
                .def("__init__", make_constructor(construct_factor_graph))
                ;
}
