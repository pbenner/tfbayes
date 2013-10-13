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

using namespace boost::python;

// library interface
// -----------------------------------------------------------------------------

char alignment_getitem(alignment_t<code_t>& a, tuple index)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= a   .length()) raise_IndexError();
        if (j >= a[i].size  ()) raise_IndexError();
        return a[i][j];
}

void alignment_setitem(alignment_t<code_t>& a, tuple index, double d)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= a   .length()) raise_IndexError();
        if (j >= a[i].size  ()) raise_IndexError();
        a[i][j] = d;
}

BOOST_PYTHON_MODULE(interface)
{
        class_<alignment_t<code_t> >("alignment_t", no_init)
                .def(init<size_t, code_t, pt_root_t>())
                .def(init<string, pt_root_t>())
                .def("__iter__",    boost::python::iterator<alignment_t<code_t> >())
                .def("__getitem__", alignment_getitem)
                .def("__setitem__", alignment_setitem)
                .def("scan",                &alignment_t<code_t>::scan<alphabet_size>)
                .def("marginal_likelihood", &alignment_t<code_t>::marginal_likelihood<alphabet_size>)
                ;
}
