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

#include <boost/noncopyable.hpp>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/operators.hpp>

#include <tfbayes/fg/distribution.hh>

using namespace boost::python;

static
double distribution_moment(const distribution_i& distribution, size_t n)
{
        switch (n) {
        case 1:
                return distribution.moment<1>();
        case 2:
                return distribution.moment<2>();
        default:
                return std::numeric_limits<double>::infinity();
        }
}

BOOST_PYTHON_MODULE(interface)
{
        class_<distribution_i, boost::noncopyable>("distribution_i", no_init)
                .def("moment", &distribution_moment)
                ;
        class_<dirac_distribution_t, bases<distribution_i> >("dirac_distribution_t", no_init)
                .def(init<double>())
                ;
        class_<exponential_family_i, bases<distribution_i>, boost::noncopyable>("exponential_family_i", no_init)
                .def("density",       &exponential_family_i::density)
                .def("base_measure",  &exponential_family_i::base_measure)
                .def("log_partition", &exponential_family_i::log_partition)
                .def("renormalize",   &exponential_family_i::renormalize)
                .def(self *= self)
                ;
        class_<normal_distribution_t, bases<exponential_family_i> >("normal_distribution_t")
                .def(init<double, double>())
                ;
        class_<gamma_distribution_t, bases<exponential_family_i> >("gamma_distribution_t")
                .def(init<double, double>())
                ;
}
