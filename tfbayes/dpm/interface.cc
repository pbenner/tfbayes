/* Copyright (C) 2011-2013 Philipp Benner
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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <tfbayes/dpm/index.hh>
#include <tfbayes/dpm/indexer.hh>
#include <tfbayes/dpm/mixture-model.hh>
#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/dpm/dpm-partition.hh>
#include <tfbayes/dpm/dpm-sampling-history.hh>

using namespace boost::python;

// tools
// -----------------------------------------------------------------------------

size_t index_i_getitem(index_i& index, size_t i)
{
        return index[i];
}

void index_i_setitem(index_i& index, size_t i, size_t d)
{
        index[i] = d;
}

template<typename T>
std::string to_string(const T& t)
{
        std::stringstream ss;
        ss << t;
        return ss.str();
}

void dpm_subset_insert(dpm_subset_t& dpm_subset, const index_i& index)
{
        dpm_subset.insert(index);
}

// interface
// -----------------------------------------------------------------------------

BOOST_PYTHON_MODULE(interface)
{
        // class definitions
        class_<sampling_history_t>("sampling_history_t")
                .def_readwrite("switches",            &sampling_history_t::switches)
                .def_readwrite("likelihood",          &sampling_history_t::likelihood)
                .def_readwrite("posterior",           &sampling_history_t::posterior)
                .def_readwrite("components",          &sampling_history_t::components)
                .def_readwrite("temperature",         &sampling_history_t::temperature)
                .def_readwrite("partitions",          &sampling_history_t::partitions)
                ;
        class_<index_i, index_i*, boost::noncopyable>("index_i", no_init)
                .def("__getitem__", index_i_getitem)
                .def("__setitem__", index_i_setitem)
                ;
        class_<index_t, bases<index_i> >("index_t")
                .def(init<size_t>())
                .def("__str__",  to_string<index_t>)
                .def("__repr__", to_string<index_t>)
                ;
        class_<seq_index_t, bases<index_i> >("seq_index_t")
                .def(init<size_t, size_t>())
                .def("__str__",  to_string<seq_index_t>)
                .def("__repr__", to_string<seq_index_t>)
                ;
        class_<dpm_subset_t>("dpm_subset_t", init<dpm_subset_tag_t>())
                .def("__iter__", boost::python::iterator<dpm_subset_t>())
                .def("__str__",  to_string<dpm_subset_t>)
                .def("__repr__", to_string<dpm_subset_t>)
                .def("__len__", &dpm_subset_t::size)
                .def("insert",  &dpm_subset_insert)
                .def("dpm_subset_tag", &dpm_subset_t::dpm_subset_tag)
                ;
        class_<dpm_partition_t>("dpm_partition_t")
                .def(vector_indexing_suite<dpm_partition_t>())
                .def("add_component", &dpm_partition_t::add_component)
                .def("__str__",  to_string<dpm_partition_t>)
                .def("__repr__", to_string<dpm_partition_t>)
                ;
        class_<dpm_partition_list_t>("dpm_partition_list_t")
                .def(vector_indexing_suite<dpm_partition_list_t>())
                .def("__str__",  to_string<dpm_partition_list_t>)
                .def("__repr__", to_string<dpm_partition_list_t>)
                ;
        class_<mixture_model_t, boost::noncopyable>("mixture_model_t", no_init)
                ;
        class_<indexer_t, boost::noncopyable>("indexer_t", no_init)
                ;
        class_<gibbs_sampler_t>("gibbs_sampler_t", no_init)
                .def(init<const mixture_model_t&,
                          const indexer_t&>())
                .def("sample", &gibbs_sampler_t::sample)
                ;
}
