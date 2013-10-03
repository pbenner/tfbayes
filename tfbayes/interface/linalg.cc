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

#include <locale>
#include <cctype>

#include <boost/python/module.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/def.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <tfbayes/utility/linalg.hh>

using namespace boost::python;

void raise_IndexError()
{
        PyErr_SetString(PyExc_IndexError, "Index out of range");
        throw error_already_set();
}

// vector class
// -----------------------------------------------------------------------------

std::vector<double>* vector_from_list(list& ns)
{
        std::vector<double>* v = new std::vector<double>(len(ns), 0);

        for (ssize_t i = 0; i < len(ns); i++) {
                (*v)[i] = extract<double>(ns[i]);
        }
        return v;
}

list vector_to_list(std::vector<double>& v)
{
        object get_iter = iterator<std::vector<double> >();
        object iter = get_iter(v);
        list l(iter);
        return l;
}

// matrix class
// -----------------------------------------------------------------------------

std::matrix<double>* matrix_from_list(list& ns)
{
        std::matrix<double>* m = new std::matrix<double>();

        for (ssize_t i = 0; i < len(ns); i++) {
                list l = extract<list>(ns[i]);
                m->push_back(std::vector<double>(len(l), 0.0));
                for (ssize_t j = 0; j < len(l); j++) {
                        (*m)[i][j] = extract<double>(l[j]);
                }
        }
        return m;
}

list matrix_to_list(std::matrix<double>& m)
{
        list l;
        for (size_t i = 0; i < m.size(); i++) {
                object get_iter = iterator<std::vector<double> >();
                object iter = get_iter(m[i]);
                l.append(list(iter));
        }
        return l;
}

double matrix_getitem(std::matrix<double>& m, tuple index)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= m   .size()) raise_IndexError();
        if (j >= m[i].size()) raise_IndexError();
        return m[i][j];
}

void matrix_setitem(std::matrix<double>& m, tuple index, double d)
{
        size_t i = extract<size_t>(index[0]);
        size_t j = extract<size_t>(index[1]);
        if (i >= m   .size()) raise_IndexError();
        if (j >= m[i].size()) raise_IndexError();
        m[i][j] = d;
}

BOOST_PYTHON_MODULE(linalg)
{
        class_<std::vector<double> >("vector")
                .def(vector_indexing_suite<std::vector<double> >() )
                .def("__init__", make_constructor(vector_from_list))
                .def(init<size_t, double>())
                .def("export", vector_to_list)
                ;

        class_<std::matrix<double> >("matrix")
                .def(vector_indexing_suite<std::matrix<double> >() )
                .def("__init__", make_constructor(matrix_from_list))
                .def("__getitem__", matrix_getitem)
                .def("__setitem__", matrix_setitem)
                .def("export", matrix_to_list)
                ;
}
