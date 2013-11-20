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
double get_moment(const exponential_family_i& distribution, size_t i)
{
        exponential_family_i::vector_t moments = distribution.moments();

        if (i >= moments.size()) {
                raise_IOError(boost::str(boost::format("Moment `%d' is not available.")
                                         % i));
        }
        return moments[i];
}

static
const exponential_family_i& access_distribution(const factor_graph_t& factor_graph, const std::string& name, size_t i)
{
        // receive distribution
        boost::optional<const exponential_family_i&> result = factor_graph.distribution(name, i);

        if (!result) {
                raise_IOError(boost::str(boost::format("Variable node `%s:%d' not found in factor graph.")
                                         % name % i));
        }
        return *result;
}
static
const exponential_family_i& access_distribution_1(const factor_graph_t& factor_graph, const std::string& name)
{
        return access_distribution(factor_graph, name, 0);
}
static
const exponential_family_i& access_distribution_2(const factor_graph_t& factor_graph, tuple t)
{
        std::string name = extract<std::string>(t[0]);
        size_t i         = extract<size_t>(t[1]);
        return access_distribution(factor_graph, name, i);
}

static
variable_node_i& access_variable_node(factor_graph_t& factor_graph, const std::string& name, size_t i)
{
        // receive distribution
        boost::optional<variable_node_i&> result = factor_graph.variable_node(name, i);

        if (!result) {
                raise_IOError(boost::str(boost::format("Variable node `%s:%d' not found in factor graph.")
                                         % name % i));
        }
        return *result;
}
static
variable_node_i& access_variable_node(factor_graph_t& factor_graph, const std::string& name)
{
        return access_variable_node(factor_graph, name, 0);
}

static
std::vector<double> call_factor_graph(factor_graph_t& factor_graph)
{
        return factor_graph();
}
static
std::vector<double> call_factor_graph(factor_graph_t& factor_graph, size_t n)
{
        return factor_graph(n);
}

// interface
// -----------------------------------------------------------------------------

BOOST_PYTHON_MODULE(interface)
{
        // distributions
        // ---------------------------------------------------------------------
        class_<exponential_family_i, boost::noncopyable>("exponential_family_i", no_init)
                .def("__call__",      &exponential_family_i::operator())
                .def("base_measure",  &exponential_family_i::base_measure)
                .def("log_partition", &exponential_family_i::log_partition)
                .def("renormalize",   &exponential_family_i::renormalize)
                .def("entropy",       &exponential_family_i::entropy)
                .def("moments",       &get_moment)
                .def(self *= self)
                ;
        class_<normal_distribution_t, bases<exponential_family_i> >("normal_distribution_t")
                .def(init<double, double>())
                ;
        class_<gamma_distribution_t, bases<exponential_family_i> >("gamma_distribution_t")
                .def(init<double, double>())
                ;
        class_<dirichlet_distribution_t, bases<exponential_family_i> >("dirichlet_distribution_t")
                .def(init<std::vector<double> >())
                ;
        class_<discrete_distribution_t, bases<exponential_family_i> >("discrete_distribution_t")
                .def(init<std::vector<double> >())
                ;
        // factor nodes
        // ---------------------------------------------------------------------
        class_<factor_node_i, boost::noncopyable>("factor_node_i", no_init)
                .def("link", static_cast<bool (factor_node_i::*)(const std::string&, variable_node_i&)>(&factor_node_i::link))
                ;
        class_<normal_fnode_t, bases<factor_node_i> >("normal_fnode_t", no_init)
                .def(init<std::string, double, double>())
                .def(init<std::string, double, double, size_t>())
                ;
        class_<gamma_fnode_t, bases<factor_node_i> >("gamma_fnode_t", no_init)
                .def(init<std::string, double, double>())
                .def(init<std::string, double, double, size_t>())
                ;
        // variable nodes
        // ---------------------------------------------------------------------
        class_<variable_node_i, boost::noncopyable>("variable_node_i", no_init)
                .def("__call__",    &variable_node_i::operator(), return_internal_reference<>())
                .def("condition",   &variable_node_i::condition)
                ;
        class_<normal_data_t, bases<variable_node_i> >("normal_data_t", no_init)
                .def(init<std::string>())
                ;
        class_<gamma_data_t, bases<variable_node_i> >("gamma_data_t", no_init)
                .def(init<std::string>())
                ;
        class_<normal_vnode_t, bases<variable_node_i> >("normal_vnode_t", no_init)
                .def(init<std::string>())
                ;
        class_<gamma_vnode_t, bases<variable_node_i> >("gamma_vnode_t", no_init)
                .def(init<std::string>())
                ;
        // factor graph
        // ---------------------------------------------------------------------
        class_<factor_graph_t>("factor_graph_t", init<>())
                .def("__call__", static_cast<std::vector<double> (*)(factor_graph_t&)>(&call_factor_graph))
                .def("__call__", static_cast<std::vector<double> (*)(factor_graph_t&, size_t)>(&call_factor_graph))
                .def("__getitem__", &access_distribution_1, return_internal_reference<>())
                .def("__getitem__", &access_distribution_2, return_internal_reference<>())
                .def("__iadd__", static_cast<factor_graph_t& (factor_graph_t::*)(const   factor_node_i&)>(&factor_graph_t::operator+=), return_internal_reference<>())
                .def("__iadd__", static_cast<factor_graph_t& (factor_graph_t::*)(const variable_node_i&)>(&factor_graph_t::operator+=), return_internal_reference<>())
                .def(self += self)
                .def("replicate", &factor_graph_t::replicate, return_internal_reference<>())
                .def("variable_node", static_cast<variable_node_i& (*)(factor_graph_t&, const std::string&)>(&access_variable_node), return_internal_reference<>())
                .def("variable_node", static_cast<variable_node_i& (*)(factor_graph_t&, const std::string&, size_t)>(&access_variable_node), return_internal_reference<>())
                .def("link", &factor_graph_t::link)
                ;
}
