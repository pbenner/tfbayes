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

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>

#include <tfbayes/phylotree/interface.hh>
#include <tfbayes/phylotree/utility.hh>
#include <tfbayes/phylotree/phylotree-parser.h>

using namespace std;

__BEGIN_DECLS

// library interface
////////////////////////////////////////////////////////////////////////////////

char* pt_print(pt_root_t* pt_root)
{
        std::ostringstream stream;
        stream << pt_root;
        std::string str = stream.str();
        const char* chr = str.c_str();

        char* result = (char*)malloc((strlen(chr)+1)*sizeof(char));
        strcpy(result, chr);

        return result;
}

const char* pt_leaf_name(pt_root_t* root, size_t leaf)
{
        return (*root)(leaf)->name.c_str();
}

size_t pt_num_leaves(pt_root_t* root)
{
        return root->n_leaves;
}

ssize_t pt_index(pt_root_t* root, const char* name)
{
        return (*root)(name)->id;
}

vector_t* pt_expectation(pt_root_t* pt_root, vector_t* observations, vector_t* prior)
{
        vector_t* result = alloc_vector(alphabet_size);
        // convert observations to an std array
        vector<code_t> tmp(observations->vec, observations->vec+observations->size);
        // compute the polynomial
        polynomial_t<code_t, alphabet_size> poly = pt_polynomial<code_t, alphabet_size>(*pt_root, tmp);

        exponent_t<code_t, alphabet_size> alpha;
        for (size_t i = 0; i < alphabet_size; i++) {
                alpha[i] = prior->vec[i];
        }

        boost::array<double, alphabet_size> expectation =
                pt_posterior_expectation<code_t, alphabet_size>(poly, alpha);

        for (size_t i = 0; i < alphabet_size; i++) {
                result->vec[i] = expectation[i];
        }

        return result;
}

vector_t* pt_approximate(pt_root_t* pt_root, vector_t* observations)
{
        vector_t* result = alloc_vector(alphabet_size);
        // convert observations to an std array
        vector<code_t> tmp(observations->vec, observations->vec+observations->size);
        // compute the polynomial
        polynomial_t<code_t, alphabet_size> poly = pt_polynomial<code_t, alphabet_size>(*pt_root, tmp);

        polynomial_t<code_t, alphabet_size> variational
                = dkl_approximate<code_t, alphabet_size>(poly);

        for (size_t i = 0; i < alphabet_size; i++) {
                result->vec[i] = variational.begin()->exponent()[i];
        }

        return result;
}

vector_t* pt_dkl_optimize(pt_root_t* pt_root, vector_t* observations)
{
        vector_t* result = alloc_vector(alphabet_size);
        // convert observations to an std array
        vector<code_t> tmp(observations->vec, observations->vec+observations->size);
        // compute the polynomial
        polynomial_t<code_t, alphabet_size> poly = pt_polynomial<code_t, alphabet_size>(*pt_root, tmp);

        exponent_t<code_t, alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;
        alpha[4] = 1;

        polynomial_t<code_t, alphabet_size> variational
                = dkl_optimize(poly, alpha);

        for (size_t i = 0; i < alphabet_size; i++) {
                result->vec[i] = variational.begin()->exponent()[i];
        }

        return result;
}

void pt_free(pt_root_t* pt_root)
{
        delete(pt_root);
}

__END_DECLS

// library interface
////////////////////////////////////////////////////////////////////////////////

#include <Python.h>

#include <locale>
#include <cctype>
#include <sstream>

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

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
                ;
        class_<pt_leaf_t, bases<pt_node_t> >("pt_leaf_t", no_init)
                ;
        class_<pt_root_t, bases<pt_node_t> >("pt_root_t", no_init)
                .def("__init__", make_constructor(parse_tree_file))
                ;
}
