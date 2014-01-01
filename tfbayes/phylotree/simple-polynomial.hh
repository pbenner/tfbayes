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

#ifndef __TFBAYES_PHYLOTREE_SIMPLE_POLYNOMIAL_HH__
#define __TFBAYES_PHYLOTREE_SIMPLE_POLYNOMIAL_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/phylotree/model.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/utility/polynomial.hh>

template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class pt_simple_polynomial_t : public polynomial_t<AS, PC> {
public:
        using polynomial_t<AS, PC>::operator=;

        pt_simple_polynomial_t(const pt_root_t& root, const std::vector<AC>& observations)
                : polynomial_t<AS, PC>() {
                operator=(likelihood_rec(root, observations)[AS]);
        }
protected:
        typedef boost::array<polynomial_t<AS, PC>, AS+1> carry_t;

        polynomial_t<AS, PC> poly_sum(const carry_t& carry) {
                polynomial_t<AS, PC> poly_sum;

                for (size_t i = 0; i < AS; i++) {
                        polynomial_term_t<AS, PC> term(1.0);
                        term.exponent()[i] = 1;
                        poly_sum += term*carry[i];
                }
                poly_sum += carry[AS];

                return poly_sum;
        }

        carry_t likelihood_leaf(const pt_leaf_t& leaf, const std::vector<AC>& observations) {
                carry_t carry;
                // i: nucleotide of ancestor
                for (size_t i = 0; i < AS; i++) {
                        carry[i] = mutation_model<AS, AC, PC>(leaf, i, observations[leaf.id]);
                }
                return carry;
        }
        carry_t likelihood_node(
                const pt_node_t& node,
                const carry_t& carry_left,
                const carry_t& carry_right,
                const std::vector<AC>& observations) {
                carry_t carry;

                // i: nucleotide of ancestor
                for (size_t i = 0; i < AS; i++) {
                        // j: nucleotide of this node
                        for (size_t j = 0; j < AS; j++) {
                                carry[i] += carry_left[j]*carry_right[j]*mutation_model<AS, AC, PC>(node, i, j);
                        }
                }
                return carry;
        }
        carry_t likelihood_root(
                const pt_node_t& node,
                const carry_t& carry_left,
                const carry_t& carry_right,
                const std::vector<AC>& observations) {
                carry_t carry;

                // j: nucleotide of this node
                for (size_t j = 0; j < AS; j++) {
                        carry[AS] += carry_left[j]*carry_right[j]*nucleotide_probability<AS, AC, PC>(j);
                }
                return carry;
        }

        carry_t likelihood_rec(const pt_node_t& node, const std::vector<AC>& observations) {
                if (node.leaf()) {
                        return likelihood_leaf(static_cast<const pt_leaf_t&>(node), observations);
                }
                else {
                        const carry_t carry_left  = likelihood_rec(node.left(),  observations);
                        const carry_t carry_right = likelihood_rec(node.right(), observations);

                        if (node.root()) {
                                return likelihood_root(node, carry_left, carry_right, observations);
                        }
                        else {
                                return likelihood_node(node, carry_left, carry_right, observations);
                        }
                }
        }
};

#endif /* __TFBAYES_PHYLOTREE_SIMPLE_POLYNOMIAL_HH__ */
