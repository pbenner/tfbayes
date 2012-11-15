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

#ifndef PHYLOTREE_APPROXIMATION_HH
#define PHYLOTREE_APPROXIMATION_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <phylotree.hh>
#include <phylotree-polynomial.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_approximate(const pt_polynomial_t<CODE_TYPE, ALPHABET_SIZE>& poly)
{
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> result;
        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> result_term(1.0);
        double norm = 0;

        /* compute normalization constant */
        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = poly.begin(); it != poly.end(); it++)
        {
                norm += it->coefficient();
        }
        /* compute approximation */
        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = poly.begin(); it != poly.end(); it++)
        {
                /* loop over the alphabet */
                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        result_term.exponent()[i] += it->coefficient()/norm * it->exponent()[i];
                }
        }
        result += result_term;

        return result;
}

#endif /* PHYLOTREE_APPROXIMATION_HH */
