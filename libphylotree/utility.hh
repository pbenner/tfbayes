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

#ifndef UTILITY_HH
#define UTILITY_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/polynomial.hh>

#include <incomplete-polynomial.hh>

#define alphabet_size 4
typedef unsigned short code_t;

std::ostream& operator<< (std::ostream& o, const exponent_t<code_t, alphabet_size>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t<code_t, alphabet_size>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t<code_t, alphabet_size>& polynomial);

size_t hash_value(const exponent_t<code_t, alphabet_size>& exponent);

#endif /* UTILITY_HH */
