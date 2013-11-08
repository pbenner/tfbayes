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

#ifndef __TFBAYES_UTILITY_DEFAULT_OPERATOR_HH__
#define __TFBAYES_UTILITY_DEFAULT_OPERATOR_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm> /* std::swap */
#include <utility>   /* std::swap */

#define default_assignment_operator(type_t) \
        type_t& operator=(const type_t& type) { \
                using std::swap; \
                type_t tmp(type); \
                swap(*this, tmp); \
                return *this; \
        }

#define virtual_assignment_operator(type_t) \
        virtual type_t& operator=(const type_t& type) { \
                using std::swap; \
                type_t tmp(type); \
                swap(*this, tmp); \
                return *this; \
        }

#define derived_assignment_operator(base_t, derived_t)  \
        virtual derived_t& operator=(const derived_t& type) { \
        std::cout << "DEFAULT ASSIGNMENT" << std::endl; \
                using std::swap; \
                derived_t tmp(static_cast<const derived_t&>(type)); \
                swap(*this, tmp); \
                return *this; \
        } \
        virtual derived_t& operator=(const base_t& type) { \
                return operator=(static_cast<const derived_t&>(type));  \
        }

#endif /* __TFBAYES_UTILITY_DEFAULT_OPERATOR_HH__ */
