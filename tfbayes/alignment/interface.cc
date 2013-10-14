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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <Python.h>

#include <locale>
#include <cctype>
#include <sstream>

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/iterator.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/interface.hh>
#include <tfbayes/phylotree/utility.hh>
#include <tfbayes/interface/exceptions.hh>
#include <tfbayes/interface/utility.hh>

using namespace boost::python;

// library interface
// -----------------------------------------------------------------------------

char alignment_getitem(alignment_t<>& a, tuple index)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= a.n_species()) raise_IndexError();
        if (j >= a.length   ()) raise_IndexError();
        return a[alignment_index_t(i,j)];
}

void alignment_setitem(alignment_t<>& a, tuple index, alphabet_code_t d)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= a.n_species()) raise_IndexError();
        if (j >= a.length   ()) raise_IndexError();
        a[alignment_index_t(i,j)] = d;
}

string alignment_str(alignment_t<>& a)
{
        return to_string(print_alignment_pretty(a));
}

BOOST_PYTHON_MODULE(interface)
{
        class_<alignment_t<> >("alignment_t", no_init)
                .def(init<size_t, pt_root_t, code_t>())
                .def(init<string, pt_root_t>())
                .def("__iter__",    boost::python::iterator<alignment_t<> >())
                .def("__getitem__", alignment_getitem)
                .def("__setitem__", alignment_setitem)
                .def("__str__",     alignment_str)
                .def("scan",                &alignment_t<>::scan<alphabet_size>)
                .def("marginal_likelihood", &alignment_t<>::marginal_likelihood<alphabet_size>)
                ;
}
