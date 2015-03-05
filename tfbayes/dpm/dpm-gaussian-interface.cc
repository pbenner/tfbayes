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

#include <tfbayes/dpm/data-gaussian.hh>
#include <tfbayes/dpm/dpm-gaussian.hh>
#include <tfbayes/dpm/dpm-partition.hh>
#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/utility/linalg.hh>

using namespace boost::python;

// tools
// -----------------------------------------------------------------------------

class dpm_gaussian_sampler_t : public gibbs_sampler_t {
public:
        dpm_gaussian_sampler_t(const dpm_gaussian_t& dpm, const data_gaussian_t& data)
                : gibbs_sampler_t(dpm, data),
                  _data          (&data)
                { }

        const data_gaussian_t& data() const {
                return *_data;
        }
        const dpm_gaussian_t& dpm() const {
                return static_cast<const dpm_gaussian_t&>(gibbs_sampler_t::dpm());
        }
private:
        const data_gaussian_t* _data;
};

size_t dpm_num_clusters(const dpm_gaussian_t& dpm) {
        return dpm.state().size();
}

dpm_partition_t dpm_partition(const dpm_gaussian_t& dpm) {
        return dpm.state().partition();
}

std::vector<double> data_get_point(const data_gaussian_t& data, const index_t& index) {
        return data[index];
}

// interface
// -----------------------------------------------------------------------------

BOOST_PYTHON_MODULE(dpm_gaussian_interface)
{
        class_<data_gaussian_t, bases<indexer_t> >("data_gaussian_t", no_init)
                .def(init<size_t,
                          const std::matrix<double>&,
                          const std::vector<double>&>())
                .def("__getitem__",                 &data_get_point)
                .def("__len__",                     &data_gaussian_t::size)
                .def("initial_means",               &data_gaussian_t::initial_means,
                     return_internal_reference<>())
                .def("initial_cluster_assignments", &data_gaussian_t::initial_cluster_assignments,
                     return_internal_reference<>())
                ;
        class_<dpm_gaussian_t, bases<mixture_model_t> >("dpm_gaussian_t", no_init)
                .def(init<double,
                          const std::matrix<double>&,
                          const std::matrix<double>&,
                          const std::vector<double>&,
                          const data_gaussian_t&>())
                .def("num_clusters", &dpm_num_clusters)
                .def("partition",    &dpm_partition)
                .def("means",        &dpm_gaussian_t::means)
                ;
        class_<dpm_gaussian_sampler_t, bases<gibbs_sampler_t> >("dpm_gaussian_sampler_t", no_init)
                .def(init<const dpm_gaussian_t&,
                          const data_gaussian_t&>())
                .def("dpm",  &dpm_gaussian_sampler_t::dpm,  return_internal_reference<>())
                .def("data", &dpm_gaussian_sampler_t::data, return_internal_reference<>())
                ;
}
