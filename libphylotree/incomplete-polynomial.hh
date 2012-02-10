/* Copyright (C) 2012 Philipp Benner
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

#ifndef INCOMPLETE_POLYNOMIAL_H
#define INCOMPLETE_POLYNOMIAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/unordered_map.hpp>

#include <phylotree.hh>
#include <utility.hh>

class node_set_t : public std::vector<pt_node_t<code_t, alphabet_size>*> {
};

class incomplete_term_t {
public:
        std::vector<node_set_t> phi;
        node_set_t incomplete;
};

class incomplete_polynomial_t : public boost::unordered_map<incomplete_term_t, double> {
};

#endif /* INCOMPLETE_POLYNOMIAL_H */
