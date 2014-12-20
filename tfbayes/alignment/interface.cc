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

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/make_constructor.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/approximation.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/phylotree/utility.hh>
#include <tfbayes/interface/exceptions.hh>
#include <tfbayes/interface/utility.hh>

using namespace boost::python;

// construct an alignment from an alignio python object
// -----------------------------------------------------------------------------

alignment_t<>* alignment_from_alignio(object a, const pt_root_t& tree)
{
        std::string type = extract<std::string>(str(a.attr("__class__")));

        /* check type of python object */
        if (type != "<class 'Bio.Align.MultipleSeqAlignment'>") {
                raise_TypeError("First argument is expected to be of type 'Bio.Align.MultipleSeqAlignment'.");
        }

        /* retrieve number of species in the alignment from the tree */
        size_t n = tree.n_leaves;

        /* start with an empty matrix */
        std::matrix<alphabet_code_t> tmp(n, 0);

        /* transfer alignment data */
        for (ssize_t i = 0; i < len(a); i++) {
                sequence_t<> s(extract<std::string>(str(a[i].attr("seq"))));
                std::string descr = extract<std::string>(str(a[i].attr("name")));
                /* the description shoud never be empty */
                assert(descr != "");
                /* extract the taxon from name */
                std::string taxon = token(descr, '.')[0];
                /* if taxon is an element of the phylogenetic tree,
                 * insert it into the alignment */
                pt_node_t::id_t id = tree.get_leaf_id(taxon);
                if (id == -1) {
                        std::cerr << boost::format("Warning: taxon `%s' not found in the phylogenetic tree.") % taxon
                                  << std::endl;
                        continue;
                }
                tmp[id] = s;
        }

        return new alignment_t<>(tmp, tree);
}

alignment_set_t<>* alignment_set_constructor(std::string filename)
{
        return new alignment_set_t<>(filename);
}

// library interface
// -----------------------------------------------------------------------------

sequence_t<> alignment_getsequence(alignment_t<>& a, std::string taxon)
{
        return a[taxon];
}

char alignment_getitem(alignment_t<>& a, boost::python::tuple index)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= a.n_species()) raise_IndexError();
        if (j >= a.length   ()) raise_IndexError();
        if (a[seq_index_t(i,j)] == a.alphabet().code('N')) {
                return 'N';
        }
        else {
                return a.alphabet().decode(a[seq_index_t(i,j)]);
        }
}

void alignment_setitem(alignment_t<>& a, boost::python::tuple index, alphabet_code_t d)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= a.n_species()) raise_IndexError();
        if (j >= a.length   ()) raise_IndexError();
        a[seq_index_t(i,j)] = a.alphabet().code(d);
}

std::vector<double> alignment_getcolumn(alignment_t<>& a, size_t i)
{
        const std::vector<alphabet_code_t> tmp(a[i]);

        return std::vector<double>(tmp.begin(), tmp.end());
}

alignment_t<> alignment_getslice(alignment_t<>& a, range_t& range)
{
        return a[range];
}

alignment_t<> alignment_set_getitem(alignment_set_t<>& a, size_t i)
{
        return a[i];
}

alignment_t<> alignment_set_getslice(alignment_set_t<>& a, range_t& range)
{
        return a[range];
}

std::string sequence_str(sequence_t<>& s)
{
        return to_string(s);
}

std::string alignment_str(alignment_t<>& a)
{
        return to_string(print_alignment_fasta(a));
}

std::string pretty_print_alignment(alignment_t<>& a)
{
        return to_string(print_alignment_pretty(a));
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
                polynomial_t<AS, PC> poly = pt_likelihood<AS, AC, PC>(tree, alignment[i]);
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
        const std::matrix<PC>& prior)
{
        std::vector<double> result;

        /* go through the alignment and compute the marginal
         * likelihood for each position */
        size_t i = 0;
        for (typename alignment_t<AC>::const_iterator it = alignment.begin();
             it != alignment.end(); it++, i++) {
                exponent_t<AS, PC> alpha(prior[i%prior.size()].begin(),
                                         prior[i%prior.size()].end  ());

                result.push_back(pt_marginal_likelihood<AS, AC, PC>
                                 (tree, *it, alpha));
        }
        return result;
}

template<typename AC, typename PC>
std::vector<double> marginal_likelihood(
        const pt_root_t& tree,
        const alignment_t<AC>& alignment,
        const std::matrix<PC>& prior)
{
        switch (alignment.alphabet().size()) {
        case 5: return marginal_likelihood<5, AC, PC>(tree, alignment, prior); break;
        default:
                std::cerr << "scan(): Invalid alphabet size."
                          << std::endl;
                exit(EXIT_FAILURE);
        }
}

template<size_t AS, typename AC, typename PC>
std::matrix<double> expectation(
        const pt_root_t& tree,
        const alignment_t<AC>& alignment,
        const std::matrix<PC>& prior)
{
        std::matrix<double> result;
        polynomial_term_t<AS> term(1.0);

        size_t i = 0;
        for (typename alignment_t<AC>::const_iterator it = alignment.begin();
             it != alignment.end(); it++, i++) {
                std::vector<double> tmp(AS, 0);
                exponent_t<AS, PC> alpha(prior[i%prior.size()].begin(),
                                         prior[i%prior.size()].end  ());

                double ml = pt_marginal_likelihood<AS, AC, PC>(tree, *it, alpha);
                pt_polynomial_t<AS, PC> poly = pt_likelihood<AS, AC, PC>(tree, *it);

                for (size_t i = 0; i < AS; i++) {
                        term.exponent()[i] = 1;
                        tmp[i] = exp(pt_marginal_likelihood(term*poly, alpha) - ml);
                        term.exponent()[i] = 0;
                }
                result.push_back(tmp);
        }

        return result;
}

template<typename AC, typename PC>
std::matrix<double> expectation(
        const pt_root_t& tree,
        const alignment_t<AC>& alignment,
        const std::matrix<PC>& prior)
{
        switch (alignment.alphabet().size()) {
        case 5: return expectation<5, AC, PC>(tree, alignment, prior); break;
        default:
                std::cerr << "scan(): Invalid alphabet size."
                          << std::endl;
                exit(EXIT_FAILURE);
        }
}

double tree_prior(
        const pt_root_t& tree,
        const boost::math::gamma_distribution<>& gamma_distribution)
{
        double result = 0.0;

        for (pt_node_t::nodes_t::const_iterator it = tree.begin_nodes();
             it != tree.end_nodes(); it++) {
                // skip the root
                if ((*it)->root()) {
                        continue;
                }
                // the gamma distribution has positive support!
                if ((*it)->d > 0) {
                        result += std::log(boost::math::pdf(gamma_distribution, (*it)->d));
                }
        }
        return result;
}

// interface
// -----------------------------------------------------------------------------

BOOST_PYTHON_MODULE(interface)
{
        class_<sequence_t<> >("sequence_t", no_init)
                .def("__str__", sequence_str)
                ;
        class_<boost::math::gamma_distribution<> >("gamma_distribution_t", no_init)
                .def(init<double, double>())
                ;
        class_<alignment_t<> >("alignment_t", no_init)
                /* the ordering of constructors is important, since
                 * alignment_from_alignio() is very general */
                .def("__init__",             make_constructor(alignment_from_alignio))
                /* this constructor needs to come last */
                .def(init<std::string, pt_root_t>())
                .def("__iter__",             boost::python::iterator<alignment_t<> >())
                .def("__getitem__",          alignment_getsequence)
                .def("__getitem__",          alignment_getitem)
                .def("__getitem__",          alignment_getslice)
                .def("__getitem__",          alignment_getcolumn)
                .def("__setitem__",          alignment_setitem)
                .def("__str__",              alignment_str)
                .def("shuffle",              &alignment_t<>::shuffle)
                .def(self += self)
                ;
        class_<alignment_set_t<> >("alignment_set_t", no_init)
                .def("__init__",             make_constructor(alignment_set_constructor))
                .def("__getitem__",          alignment_set_getitem)
                .def("__getitem__",          alignment_set_getslice)
                .def("__len__",              &alignment_set_t<>::size)
                ;
                
        // functions on alignments
        def("approximate",            &approximate        <alphabet_code_t, double>);
        def("marginal_likelihood",    &marginal_likelihood<alphabet_code_t, double>);
        def("scan",                   &scan               <alphabet_code_t, double>);
        def("expectation",            &expectation        <alphabet_code_t, double>);
        def("tree_prior",             &tree_prior                                  );
        def("pretty_print_alignment", &pretty_print_alignment                      );
}
