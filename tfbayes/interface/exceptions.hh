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


#ifndef __TFBAYES_INTERFACE_INTERFACE_EXCEPTIONS_HH__
#define __TFBAYES_INTERFACE_INTERFACE_EXCEPTIONS_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <Python.h>

#include <string>
#include <boost/format.hpp>
#include <boost/python/errors.hpp>

static inline
void raise_IndexError()
{
        PyErr_SetString(PyExc_IndexError, "Index out of range");
        throw boost::python::error_already_set();
}

static inline
void raise_IOError(const std::string& msg)
{
        PyErr_SetString(PyExc_IOError, msg.c_str());
        throw boost::python::error_already_set();
}

static inline
void raise_IOError(boost::basic_format<char>& msg)
{
        raise_IOError(boost::str(msg));
}

static inline
void raise_TypeError(const std::string& msg)
{
        PyErr_SetString(PyExc_TypeError, msg.c_str());
        throw boost::python::error_already_set();
}

#endif /* __TFBAYES_INTERFACE_INTERFACE_EXCEPTIONS_HH__ */
