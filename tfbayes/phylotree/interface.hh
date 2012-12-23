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

#ifndef INTERFACE_HH
#define INTERFACE_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/phylotree/phylotree-approximation.hh>
#include <tfbayes/phylotree/posterior.hh>
#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/exception/exception.h>
#include <tfbayes/utility/linalg.h>

#define alphabet_size 5
//typedef unsigned short code_t;
typedef float code_t;

#define EXTERN_C_START extern "C" {
#define EXTERN_C_END   }

#endif /* INTERFACE_HH */
