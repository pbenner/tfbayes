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
#include <tfbayes/dpm/dpm-partition.hh>
#include <tfbayes/dpm/dpm-sampling-history.hh>
#include <tfbayes/dpm/dpm-tfbs-options.hh>
#include <tfbayes/dpm/dpm-tfbs-sampler.hh>

using namespace boost::python;

// tools
// -----------------------------------------------------------------------------

void baseline_priors_push_back(
        baseline_priors_t& baseline_priors,
        const std::matrix<double>& m)
{
        baseline_priors.push_back(m);
}

void baseline_names_push_back(
        baseline_names_t& baseline_names,
        const std::string& s)
{
        baseline_names.push_back(s);
}

// interface
// -----------------------------------------------------------------------------

BOOST_PYTHON_MODULE(dpm_tfbs_interface)
{
        class_<tfbs_options_t>("tfbs_options_t")
                .def_readwrite("phylogenetic_file",    &tfbs_options_t::phylogenetic_file)
                .def_readwrite("alignment_file",       &tfbs_options_t::alignment_file)
                .def_readwrite("alpha",                &tfbs_options_t::alpha)
                .def_readwrite("discount",             &tfbs_options_t::discount)
                .def_readwrite("_lambda_",             &tfbs_options_t::lambda)
                .def_readwrite("block_samples",        &tfbs_options_t::block_samples)
                .def_readwrite("block_samples_period", &tfbs_options_t::block_samples_period)
                .def_readwrite("metropolis_proposals", &tfbs_options_t::metropolis_proposals)
                .def_readwrite("optimize",             &tfbs_options_t::optimize)
                .def_readwrite("optimize_period",      &tfbs_options_t::optimize_period)
                .def_readwrite("initial_temperature",  &tfbs_options_t::initial_temperature)
                .def_readwrite("process_prior",        &tfbs_options_t::process_prior)
                .def_readwrite("background_model",     &tfbs_options_t::background_model)
                .def_readwrite("background_alpha",     &tfbs_options_t::background_alpha)
                .def_readwrite("background_context",   &tfbs_options_t::background_context)
                .def_readwrite("background_beta",      &tfbs_options_t::background_beta)
                .def_readwrite("background_gamma",     &tfbs_options_t::background_gamma)
                .def_readwrite("background_cache",     &tfbs_options_t::background_cache)
                .def_readwrite("background_weights",   &tfbs_options_t::background_weights)
                .def_readwrite("baseline_names",       &tfbs_options_t::baseline_names)
                .def_readwrite("baseline_lengths",     &tfbs_options_t::baseline_lengths)
                .def_readwrite("baseline_priors",      &tfbs_options_t::baseline_priors)
                .def_readwrite("baseline_weights",     &tfbs_options_t::baseline_weights)
                .def_readwrite("population_size",      &tfbs_options_t::population_size)
                .def_readwrite("threads",              &tfbs_options_t::threads)
                .def_readwrite("socket_file",          &tfbs_options_t::socket_file)
                .def_readwrite("verbose",              &tfbs_options_t::verbose)
                ;
        class_<baseline_names_t>("baseline_names_t")
                .def("__iter__", boost::python::iterator<baseline_names_t>())
                .def("append",  &baseline_names_push_back)
                ;
        class_<baseline_priors_t>("baseline_priors_t")
                .def("__iter__", boost::python::iterator<baseline_priors_t>())
                .def("append", &baseline_priors_push_back)
                ;
        class_<data_tfbs_t>("data_tfbs_t", no_init)
                .def(init<const std::string&>())
                ;
        class_<dpm_tfbs_t>("dpm_tfbs_t", no_init)
                .def(init<const tfbs_options_t&, const data_tfbs_t&>())
                .def("map",     &dpm_tfbs_t::map)
                .def("mean",    &dpm_tfbs_t::mean)
                .def("median",  &dpm_tfbs_t::median)
                ;
        class_<dpm_tfbs_pmcmc_t, boost::noncopyable>("dpm_tfbs_pmcmc_t", no_init)
                .def(init<tfbs_options_t>())
                .def(init<tfbs_options_t, sampling_history_t>())
                .def("__call__", &dpm_tfbs_pmcmc_t::operator())
                .def("save",     &dpm_tfbs_pmcmc_t::save)
                ;
}
