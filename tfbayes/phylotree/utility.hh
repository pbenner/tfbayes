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

std::ostream& operator<< (std::ostream& o, const exponent_t        <4, char>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t <4, char>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t      <4, char>& polynomial);
std::ostream& operator<< (std::ostream& o, const exponent_t        <5, char>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t <5, char>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t      <5, char>& polynomial);

std::ostream& operator<< (std::ostream& o, const exponent_t        <4, short>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t <4, short>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t      <4, short>& polynomial);
std::ostream& operator<< (std::ostream& o, const exponent_t        <5, short>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t <5, short>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t      <5, short>& polynomial);

std::ostream& operator<< (std::ostream& o, const exponent_t        <4, float>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t <4, float>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t      <4, float>& polynomial);
std::ostream& operator<< (std::ostream& o, const exponent_t        <5, float>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t <5, float>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t      <5, float>& polynomial);

std::ostream& operator<< (std::ostream& o, const exponent_t        <4, double>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t <4, double>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t      <4, double>& polynomial);
std::ostream& operator<< (std::ostream& o, const exponent_t        <5, double>& exponent);
std::ostream& operator<< (std::ostream& o, const polynomial_term_t <5, double>& term);
std::ostream& operator<< (std::ostream& o, const polynomial_t      <5, double>& polynomial);

size_t hash_value(const exponent_t< 2, char>& exponent);
size_t hash_value(const exponent_t< 3, char>& exponent);
size_t hash_value(const exponent_t< 4, char>& exponent);
size_t hash_value(const exponent_t< 5, char>& exponent);
size_t hash_value(const exponent_t<10, char>& exponent);

size_t hash_value(const exponent_t< 2, short>& exponent);
size_t hash_value(const exponent_t< 3, short>& exponent);
size_t hash_value(const exponent_t< 4, short>& exponent);
size_t hash_value(const exponent_t< 5, short>& exponent);
size_t hash_value(const exponent_t<10, short>& exponent);

size_t hash_value(const exponent_t< 2, float>& exponent);
size_t hash_value(const exponent_t< 3, float>& exponent);
size_t hash_value(const exponent_t< 4, float>& exponent);
size_t hash_value(const exponent_t< 5, float>& exponent);
size_t hash_value(const exponent_t<10, float>& exponent);

size_t hash_value(const exponent_t< 2, double>& exponent);
size_t hash_value(const exponent_t< 3, double>& exponent);
size_t hash_value(const exponent_t< 4, double>& exponent);
size_t hash_value(const exponent_t< 5, double>& exponent);
size_t hash_value(const exponent_t<10, double>& exponent);

#endif /* UTILITY_HH */
