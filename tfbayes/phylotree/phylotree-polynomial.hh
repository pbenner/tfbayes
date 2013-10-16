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

#ifndef PHYLOTREE_POLYNOMIAL_HH
#define PHYLOTREE_POLYNOMIAL_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/polynomial.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

template <size_t AS, typename AC = alphabet_code_t, typename PC = short>
class pt_polynomial_t : public polynomial_t<AS, PC> {
public:
        using polynomial_t<AS, PC>::operator=;

        pt_polynomial_t(const pt_root_t& root, const std::vector<AC>& observations) {
                if (root.outgroup() && observations[root.outgroup()->id] != -1) {
                        operator=(poly_sum(likelihood_rec(root, observations), *root.outgroup(), observations));
                }
                else {
                        operator=(poly_sum(likelihood_rec(root, observations)));
                }
        }

private:
        /******************************************************************************
         * typedefs
         ******************************************************************************/

        typedef boost::array<polynomial_t<AS, PC>, AS+1> carry_t;

        /******************************************************************************
         * Algorithm by PY and PB
         ******************************************************************************/

        static polynomial_t<AS, PC> poly_sum(const carry_t& carry) {
                polynomial_t<AS, PC> poly_sum;

                for (size_t i = 0; i < AS; i++) {
                        poly_sum += nucleotide_probability(i)*carry[i];
                }
                poly_sum += carry[AS];

                return poly_sum;
        }
        static polynomial_t<AS, PC> poly_sum(const carry_t& carry, const pt_leaf_t& outgroup,
                                                               const std::vector<AC>& observations) {
                polynomial_t<AS, PC> poly_sum;

                for (size_t i = 0; i < AS; i++) {
                        poly_sum += mutation_model(outgroup, observations[outgroup.id], i)*carry[i];
                }
                poly_sum += carry[AS];
                poly_sum *= nucleotide_probability(observations[outgroup.id]);

                return poly_sum;
        }
        static carry_t likelihood_leaf(const pt_node_t& node, const std::vector<AC>& observations) {
                carry_t carry;
                const pt_leaf_t& leaf = static_cast<const pt_leaf_t&>(node);
                if (observations[leaf.id] == -1) {
                        // this leaf should be ignored, no nucleotide
                        // is present (missing data)
                        carry[AS] = 1.0;
                }
                else {
                        // leaf is present and has a nucleotide
                        carry[observations[leaf.id]] += 1.0;
                }
                return carry;
        }
        static carry_t likelihood_node(
                const pt_node_t& node,
                const carry_t& carry_left,
                const carry_t& carry_right,
                const std::vector<AC>& observations) {
                double pm_left  = 1.0-exp(-node.left ().d);
                double pm_right = 1.0-exp(-node.right().d);
                carry_t carry;
                const polynomial_t<AS, PC> poly_sum_left  = poly_sum(carry_left);
                const polynomial_t<AS, PC> poly_sum_right = poly_sum(carry_right);

                carry[AS] +=
                        ((1.0-pm_left )*carry_left [AS] + pm_left *poly_sum_left)*
                        ((1.0-pm_right)*carry_right[AS] + pm_right*poly_sum_right);

                for (size_t i = 0; i < AS; i++) {
                        carry[i] += (1.0-pm_left )*(1.0-pm_right)*carry_left [i]*carry_right[AS];
                        carry[i] += (1.0-pm_left )*     pm_right *carry_left [i]*poly_sum_right;
                        carry[i] += (1.0-pm_right)*(1.0-pm_left )*carry_right[i]*carry_left [AS];
                        carry[i] += (1.0-pm_right)*     pm_left  *carry_right[i]*poly_sum_left;
                        carry[i] += (1.0-pm_left )*(1.0-pm_right)*carry_left [i]*carry_right[i];
                }
                return carry;
        }
        static carry_t likelihood_rec(const pt_node_t& node, const std::vector<AC>& observations) {
                if (node.leaf()) {
                        return likelihood_leaf(node, observations);
                }
                else {
                        const carry_t carry_left  = likelihood_rec(node.left (),  observations);
                        const carry_t carry_right = likelihood_rec(node.right(), observations);

                        return likelihood_node(node, carry_left, carry_right, observations);
                }
        }
        static polynomial_term_t<AS, PC> nucleotide_probability(AC x) {
                polynomial_term_t<AS, PC> px(1.0);
                px.exponent()[x] = 1;
                return px;
        }
        static polynomial_t<AS, PC> mutation_model(const pt_node_t& node, AC x, AC y) {
                polynomial_t<AS, PC> poly;
                polynomial_term_t<AS, PC> py = nucleotide_probability(y);
                poly += node.mutation_probability()*py;
                if (x == y) {
                        poly += (1.0-node.mutation_probability());
                }
                return poly;
        }
};

#endif /* PHYLOTREE_POLYNOMIAL_HH */
