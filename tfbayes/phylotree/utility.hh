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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/utility/polynomial.hh>

#include <tfbayes/phylotree/incomplete-expression.hh>

std::ostream& operator<< (std::ostream& o, const exponent_t<short, 4>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t<short, 4>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t<short, 4>& polynomial);
std::ostream& operator<< (std::ostream& o, const exponent_t<short, 5>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t<short, 5>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t<short, 5>& polynomial);

std::ostream& operator<< (std::ostream& o, const exponent_t<float, 4>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t<float, 4>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t<float, 4>& polynomial);
std::ostream& operator<< (std::ostream& o, const exponent_t<float, 5>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t<float, 5>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t<float, 5>& polynomial);

std::ostream& operator<< (std::ostream& o, const exponent_t<double, 4>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t<double, 4>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t<double, 4>& polynomial);
std::ostream& operator<< (std::ostream& o, const exponent_t<double, 5>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t<double, 5>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t<double, 5>& polynomial);

size_t hash_value(const exponent_t<short, 2>& exponent);
size_t hash_value(const exponent_t<short, 3>& exponent);
size_t hash_value(const exponent_t<short, 4>& exponent);
size_t hash_value(const exponent_t<short, 5>& exponent);
size_t hash_value(const exponent_t<short, 10>& exponent);

size_t hash_value(const exponent_t<float, 2>& exponent);
size_t hash_value(const exponent_t<float, 3>& exponent);
size_t hash_value(const exponent_t<float, 4>& exponent);
size_t hash_value(const exponent_t<float, 5>& exponent);
size_t hash_value(const exponent_t<float, 10>& exponent);

size_t hash_value(const exponent_t<double, 2>& exponent);
size_t hash_value(const exponent_t<double, 3>& exponent);
size_t hash_value(const exponent_t<double, 4>& exponent);
size_t hash_value(const exponent_t<double, 5>& exponent);
size_t hash_value(const exponent_t<double, 10>& exponent);

#endif /* UTILITY_HH */
