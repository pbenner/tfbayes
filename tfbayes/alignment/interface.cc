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

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/make_constructor.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/interface.hh>
#include <tfbayes/phylotree/utility.hh>
#include <tfbayes/interface/exceptions.hh>
#include <tfbayes/interface/utility.hh>

using namespace boost::python;

// construct an alignment from an alignio python object
// -----------------------------------------------------------------------------

#include <string>

alignment_t<>* alignment_from_alignio(object a, const pt_root_t& tree)
{
        string type = extract<string>(str(a.attr("__class__")));

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
                sequence_t<> s(extract<string>(str(a[i].attr("seq"))));
                string descr = extract<string>(str(a[i].attr("name")));
                /* the description shoud never be empty */
                assert(descr != "");
                /* extract the taxon from name */
                string taxon = token(descr, '.')[0];
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

// library interface
// -----------------------------------------------------------------------------

sequence_t<> alignment_getsequence(alignment_t<>& a, string taxon)
{
        return a[taxon];
}

char alignment_getitem(alignment_t<>& a, tuple index)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= a.n_species()) raise_IndexError();
        if (j >= a.length   ()) raise_IndexError();
        if (a[alignment_index_t(i,j)] == a.alphabet().code('N')) {
                return 'N';
        }
        else {
                return a.alphabet().decode(a[alignment_index_t(i,j)]);
        }
}

void alignment_setitem(alignment_t<>& a, tuple index, alphabet_code_t d)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= a.n_species()) raise_IndexError();
        if (j >= a.length   ()) raise_IndexError();
        a[alignment_index_t(i,j)] = a.alphabet().code(d);
}

string sequence_str(sequence_t<>& s)
{
        return to_string(s);
}

string alignment_str(alignment_t<>& a)
{
        return to_string(print_alignment_pretty(a));
}

BOOST_PYTHON_MODULE(interface)
{
        class_<sequence_t<> >("sequence_t", no_init)
                .def("__str__", sequence_str)
                ;
        class_<alignment_t<> >("alignment_t", no_init)
                .def(init<string, pt_root_t>())
                .def("__init__",    make_constructor(alignment_from_alignio))
                .def("__iter__",    boost::python::iterator<alignment_t<> >())
                .def("__getitem__", alignment_getsequence)
                .def("__getitem__", alignment_getitem)
                .def("__setitem__", alignment_setitem)
                .def("__str__",     alignment_str)
                .def("scan",                &alignment_t<>::scan<alphabet_size>)
                .def("marginal_likelihood", &alignment_t<>::marginal_likelihood<alphabet_size>)
                ;
}
