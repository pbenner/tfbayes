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

#include <tfbayes/dpm/init.hh>
#include <tfbayes/dpm/index.hh>
#include <tfbayes/dpm/dpm-partition.hh>
#include <tfbayes/dpm/dpm-sampling-history.hh>
#include <tfbayes/dpm/dpm-tfbs-options.hh>
#include <tfbayes/dpm/dpm-tfbs-sampler.hh>

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

void baseline_priors_push_back(
        baseline_priors_t& baseline_priors,
        const std::matrix<double>& m)
{
        baseline_priors.push_back(m);
}

void baseline_tags_push_back(
        baseline_tags_t& baseline_tags,
        const std::string& s)
{
        baseline_tags.push_back(s);
}

template<typename T>
std::string to_string(const T& t)
{
        std::stringstream ss;
        ss << t;
        return ss.str();
}

// interface
// -----------------------------------------------------------------------------

BOOST_PYTHON_MODULE(dpm_tfbs_interface)
{
        // init dpm library
        def("dpm_init", __dpm_init__);
        def("dpm_free", __dpm_free__);
        // class definitions
        class_<sampling_history_t>("sampling_history_t")
                .def_readwrite("switches",            &sampling_history_t::switches)
                .def_readwrite("likelihood",          &sampling_history_t::likelihood)
                .def_readwrite("posterior",           &sampling_history_t::posterior)
                .def_readwrite("components",          &sampling_history_t::components)
                .def_readwrite("temperature",         &sampling_history_t::temperature)
                .def_readwrite("partitions",          &sampling_history_t::partitions)
                ;
        class_<tfbs_options_t>("tfbs_options_t")
                .def_readwrite("phylogenetic_file",   &tfbs_options_t::phylogenetic_file)
                .def_readwrite("alignment_file",      &tfbs_options_t::alignment_file)
                .def_readwrite("tfbs_length",         &tfbs_options_t::tfbs_length)
                .def_readwrite("alpha",               &tfbs_options_t::alpha)
                .def_readwrite("discount",            &tfbs_options_t::discount)
                .def_readwrite("_lambda_",            &tfbs_options_t::lambda)
                .def_readwrite("initial_temperature", &tfbs_options_t::initial_temperature)
                .def_readwrite("process_prior",       &tfbs_options_t::process_prior)
                .def_readwrite("background_model",    &tfbs_options_t::background_model)
                .def_readwrite("background_alpha",    &tfbs_options_t::background_alpha)
                .def_readwrite("background_context",  &tfbs_options_t::background_context)
                .def_readwrite("background_weights",  &tfbs_options_t::background_weights)
                .def_readwrite("baseline_priors",     &tfbs_options_t::baseline_priors)
                .def_readwrite("baseline_weights",    &tfbs_options_t::baseline_weights)
                .def_readwrite("baseline_tags",       &tfbs_options_t::baseline_tags)
                .def_readwrite("population_size",     &tfbs_options_t::population_size)
                .def_readwrite("socket_file",         &tfbs_options_t::socket_file)
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
                .def("insert", &dpm_subset_t::insert)
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
        class_<baseline_priors_t>("baseline_priors_t")
                .def("__iter__", boost::python::iterator<baseline_priors_t>())
                .def("append", &baseline_priors_push_back)
                ;
        class_<baseline_tags_t>("baseline_tags_t")
                .def("__iter__", boost::python::iterator<baseline_tags_t>())
                .def("append",  &baseline_tags_push_back)
                ;
        class_<dpm_tfbs_pmcmc_t, boost::noncopyable>("dpm_tfbs_pmcmc_t", no_init)
                .def(init<tfbs_options_t>())
                .def(init<tfbs_options_t, sampling_history_t>())
                .def("sample", &dpm_tfbs_pmcmc_t::sample)
                .def("save",   &dpm_tfbs_pmcmc_t::save)
                ;
}
