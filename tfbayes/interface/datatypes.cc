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
#include <list>
#include <string>

#include <boost/python/module.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/init.hpp>
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

// implicit type conversion
// -----------------------------------------------------------------------------

struct iterable_converter
{
        /// @note Registers converter from a python interable type to the
        ///       provided type.
        template <typename Container>
        iterable_converter&
        from_python() {
                boost::python::converter::registry::push_back(
                        &iterable_converter::convertible,
                        &iterable_converter::construct<Container>,
                        boost::python::type_id<Container>());
                return *this;
        }
        static void* convertible(PyObject* object) {
                return PyObject_GetIter(object) ? object : NULL;
        }
        template <typename Container>
        static void construct(
                PyObject* object,
                boost::python::converter::rvalue_from_python_stage1_data* data) {
                namespace python = boost::python;
                // Object is a borrowed reference, so create a handle indicting it is
                // borrowed for proper reference counting.
                python::handle<> handle(python::borrowed(object));
                
                // Obtain a handle to the memory block that the converter has allocated
                // for the C++ type.
                typedef python::converter::rvalue_from_python_storage<Container>
                        storage_type;
                void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;
                
                typedef python::stl_input_iterator<typename Container::value_type>
                        iterator;
                
                // Allocate the C++ type into the converter's memory block, and assign
                // its handle to the converter's convertible variable.  The C++
                // container is populated by passing the begin and end iterators of
                // the python object to the container's constructor.
                data->convertible = new (storage) Container(
                        iterator(python::object(handle)), // begin
                        iterator());                      // end
        }
};

// register classes
// -----------------------------------------------------------------------------

BOOST_PYTHON_MODULE(datatypes)
{
        class_<std::vector<double> >("vector")
                .def(vector_indexing_suite<std::vector<double> >() )
                .def("__init__", make_constructor(vector_from_list))
                .def("export", vector_to_list)
                ;
        class_<std::matrix<double> >("matrix")
                .def(vector_indexing_suite<std::matrix<double> >() )
                .def("__init__", make_constructor(matrix_from_list))
                .def("__getitem__", matrix_getitem)
                .def("__setitem__", matrix_setitem)
                .def("export", matrix_to_list)
                ;
        iterable_converter()
              .from_python<std::vector<double> >()
              .from_python<std::matrix<double> >()
              .from_python<std::list<std::string> >()
              .from_python<std::list<std::matrix<double > > >()
              ;
}
