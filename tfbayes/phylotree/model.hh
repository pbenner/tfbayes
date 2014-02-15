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

#ifndef __TFBAYES_PHYLOTREE_MODEL_HH__
#define __TFBAYES_PHYLOTREE_MODEL_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/utility/polynomial.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

template <size_t AS, typename AC, typename PC>
polynomial_term_t<AS, PC> nucleotide_probability(AC x) {
        polynomial_term_t<AS, PC> px(1.0);
        px.exponent()[x] = 1;
        return px;
}

template <size_t AS, typename AC, typename PC>
polynomial_t<AS, PC> mutation_model(const pt_node_t& node, AC x, AC y) {
        polynomial_t<AS, PC> poly;
        polynomial_term_t<AS, PC> py = nucleotide_probability<AS, AC, PC>(y);
        poly += node.mutation_probability()*py;
        if (x == y) {
                poly += (1.0-node.mutation_probability());
        }
        return poly;
}

#endif /* __TFBAYES_PHYLOTREE_MODEL_HH__ */
