/* Copyright (C) 2012 Philipp Benner
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

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <tfbayes/phylotree/parser.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/treespace.hh>
#include <tfbayes/interface/exceptions.hh>
#include <tfbayes/interface/utility.hh>

using namespace boost::python;

// tools
// -----------------------------------------------------------------------------

std::string print_tree(const pt_root_t& tree)
{
        return to_string(newick_format(tree));
}

pt_root_t* parse_tree_file(const std::string& filename)
{
        std::list<pt_root_t> tree_list = parse_tree_list(filename);

        if (tree_list.size() != 1) {
                raise_IOError(boost::format("%s: File must contain a single phylogenetic tree.")
                              % filename);
        }
        return new pt_root_t(tree_list.front());
}

pt_root_t* parse_tree_file_ref(const std::string& filename, const pt_root_t& ref_tree)
{
        std::list<pt_root_t> tree_list = parse_tree_list(filename, 0, 1, ref_tree);

        if (tree_list.size() != 1) {
                raise_IOError(boost::format("%s: File must contain a single phylogenetic tree.")
                              % filename);
        }
        return new pt_root_t(tree_list.front());
}

pt_root_t geodesic_call(geodesic_t geodesic, double lambda)
{
        return geodesic(lambda).export_tree();
}

// interface
// -----------------------------------------------------------------------------

BOOST_PYTHON_MODULE(interface)
{
        class_<pt_node_t>("pt_node_t", no_init)
                .def("__str__", print_tree)
                .def_readonly("n_leaves", &pt_node_t::n_leaves)
                .def_readonly("n_nodes",  &pt_node_t::n_nodes)
                .def_readonly("d",        &pt_node_t::d)
                .def("leaf",              &pt_node_t::leaf)
                .def("root",              &pt_node_t::root)
                .def("scale",             &pt_node_t::scale)
                ;
        class_<pt_leaf_t, bases<pt_node_t> >("pt_leaf_t", no_init)
                ;
        class_<pt_root_t, bases<pt_node_t> >("pt_root_t", no_init)
                .def("__init__",    make_constructor(parse_tree_file))
                .def("__init__",    make_constructor(parse_tree_file_ref))
                .def("get_node_id", &pt_root_t::get_node_id)
                .def("get_leaf_id", &pt_root_t::get_leaf_id)
                ;
        class_<geodesic_t>("geodesic_t", no_init)
                .def(init<pt_root_t, pt_root_t>())
                .def("__call__", &geodesic_call)
                .def("length",   &geodesic_t::length)
                ;
}
