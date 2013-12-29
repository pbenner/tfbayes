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

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree-approximation.hh>
#include <tfbayes/phylotree/phylotree-parser.h>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
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

// alignment functions
// -----------------------------------------------------------------------------

template<size_t AS, typename AC, typename PC>
std::matrix<double> approximate(
        const pt_root_t& tree,
        const alignment_t<AC>& alignment)
{
        std::matrix<double> result(alignment.length(), AS);

        for (size_t i = 0; i < alignment.length(); i++) {
                // compute the polynomial
                polynomial_t<AS, PC> poly = pt_polynomial_t<AS, AC, PC>(tree, alignment[i]);
                polynomial_t<AS, PC> variational
                        = dkl_approximate<AS, PC>(poly);

                for (size_t j = 0; j < AS; j++) {
                        result[i][j] = variational.begin()->exponent()[j];
                }
        }
        return result;
}

template<typename AC, typename PC>
std::matrix<double> approximate(
        const pt_root_t& tree,
        const alignment_t<AC>& alignment)
{
        switch (alignment.alphabet().size()) {
        case 5: return approximate<5, AC, PC>(tree, alignment); break;
        default: 
                std::cerr << "scan(): Invalid alphabet size."
                          << std::endl;
                exit(EXIT_FAILURE);
        }
}

template<size_t AS, typename AC, typename PC>
std::vector<double> scan(
        const pt_root_t& tree,
        const alignment_t<AC>& alignment,
        std::matrix<PC>& counts)
{
        std::vector<double> result(alignment.length(), 0);
        std::vector<exponent_t<AS, PC> > exponents;

        for (size_t j = 0; j < counts.size(); j++) {
                exponent_t<AS, PC> tmp
                        (counts[j].begin(), counts[j].end());
                exponents.push_back(tmp);
        }

        for (typename alignment_t<AC>::const_iterator it = alignment.begin();
             it != alignment.end(); it++) {
                result[it - alignment.begin()] = 0;
                if (it - alignment.begin() + counts.size() > alignment.length()) {
                        // do not exit the loop here so that every
                        // position is initialized
                        continue;
                }
                for (typename alignment_t<AC>::const_iterator is(it); is < it + counts.size(); is++) {
                        result[it - alignment.begin()] += pt_marginal_likelihood<AS, AC, PC>(
                                tree, *is, exponents[is-it]);
                }
        }
        return result;
}

template<typename AC, typename PC>
std::vector<double> scan(
        const pt_root_t& tree,
        const alignment_t<AC>& alignment,
        std::matrix<PC>& counts)
{
        switch (alignment.alphabet().size()) {
        case 5: return scan<5, AC, PC>(tree, alignment, counts); break;
        default: 
                std::cerr << "scan(): Invalid alphabet size."
                          << std::endl;
                exit(EXIT_FAILURE);
        }
}

template<size_t AS, typename AC, typename PC>
std::vector<double> marginal_likelihood(
        const pt_root_t& tree,
        const alignment_t<AC>& alignment,
        const std::vector<PC>& prior)
{
        exponent_t<AS, PC> alpha(prior.begin(), prior.end());
        std::vector<double> result;

        /* go through the alignment and compute the marginal
         * likelihood for each position */
        for (typename alignment_t<AC>::const_iterator it = alignment.begin();
             it != alignment.end(); it++) {
                result.push_back(pt_marginal_likelihood<AS, AC, PC>
                                 (tree, *it, alpha));
        }
        return result;
}

template<typename AC, typename PC>
std::vector<double> marginal_likelihood(
        const pt_root_t& tree,
        const alignment_t<AC>& alignment,
        const std::vector<PC>& prior)
{
        switch (alignment.alphabet().size()) {
        case 5: return marginal_likelihood<5, AC, PC>(tree, alignment, prior); break;
        default:
                std::cerr << "scan(): Invalid alphabet size."
                          << std::endl;
                exit(EXIT_FAILURE);
        }
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
                .def("__init__",    make_constructor(parse_tree_file))
                .def("get_node_id", &pt_root_t::get_node_id)
                .def("get_leaf_id", &pt_root_t::get_leaf_id)
                ;
        // functions on alignments
        def("approximate",         &approximate        <alphabet_code_t, double>);
        def("marginal_likelihood", &marginal_likelihood<alphabet_code_t, double>);
        def("scan",                &scan               <alphabet_code_t, double>);
}
