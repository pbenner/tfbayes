/* Copyright (C) 2012-2013 Philipp Benner
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

#ifndef __TFBAYES_PHYLOTREE_POLYNOMIAL_HH__
#define __TFBAYES_PHYLOTREE_POLYNOMIAL_HH__

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

namespace tfbayes_detail {
        /******************************************************************************
         * typedefs
         ******************************************************************************/

        template <size_t AS, typename AC, typename PC>
        class carry_t : public boost::array<polynomial_t<AS, PC>, AS+1>
        { };

        /******************************************************************************
         * Algorithm by PY and PB
         ******************************************************************************/

        template <size_t AS, typename AC, typename PC>
        polynomial_term_t<AS, PC> nucleotide_probability(AC x) {
                polynomial_term_t<AS, PC> px(1.0);
                px.exponent()[x] = 1;
                return px;
        }
        template <size_t AS, typename AC, typename PC>
        polynomial_t<AS, PC> poly_sum(const carry_t<AS, AC, PC>& carry) {
                polynomial_t<AS, PC> poly_sum;

                for (size_t i = 0; i < AS; i++) {
                        poly_sum += nucleotide_probability<AS, AC, PC>(i)*carry[i];
                }
                poly_sum += carry[AS];

                return poly_sum;
        }
        template <size_t AS, typename AC, typename PC>
        carry_t<AS, AC, PC> likelihood_leaf(const pt_node_t& node, const std::vector<AC>& observations) {
                carry_t<AS, AC, PC> carry;
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
        template <size_t AS, typename AC, typename PC>
        carry_t<AS, AC, PC> likelihood_node(
                const pt_node_t& node,
                const carry_t<AS, AC, PC>& carry_left,
                const carry_t<AS, AC, PC>& carry_right,
                const std::vector<AC>& observations) {
                double pm_left  = 1.0-exp(-node.left ().d);
                double pm_right = 1.0-exp(-node.right().d);
                carry_t<AS, AC, PC> carry;
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
        template <size_t AS, typename AC, typename PC>
        carry_t<AS, AC, PC> likelihood_rec(const pt_node_t& node, const std::vector<AC>& observations) {
                if (node.leaf()) {
                        return likelihood_leaf<AS, AC, PC>(node, observations);
                }
                else {
                        const carry_t<AS, AC, PC> carry_left  = likelihood_rec<AS, AC, PC>(node.left (), observations);
                        const carry_t<AS, AC, PC> carry_right = likelihood_rec<AS, AC, PC>(node.right(), observations);

                        return likelihood_node<AS, AC, PC>(node, carry_left, carry_right, observations);
                }
        }
        template <size_t AS, typename AC, typename PC>
        carry_t<AS, AC, PC> likelihood_root(
                const pt_root_t& root,
                const std::vector<AC>& observations) {

                // if the root has no outgroup then treat it as a
                // simple node
                if (!root.outgroup()) {
                        return likelihood_rec<AS, AC, PC>(root, observations);
                }
                const carry_t<AS, AC, PC> carry_left  = likelihood_rec<AS, AC, PC>( root,            observations);
                const carry_t<AS, AC, PC> carry_right = likelihood_rec<AS, AC, PC>(*root.outgroup(), observations);

                const polynomial_t<AS, PC> poly_sum_left  = poly_sum(carry_left);
                const polynomial_t<AS, PC> poly_sum_right = poly_sum(carry_right);

                double pm_right = 1.0-exp(-root.outgroup()->d);
                carry_t<AS, AC, PC> carry;

                carry[AS] +=
                        carry_left [AS]*((1.0-pm_right)*carry_right[AS] + pm_right*poly_sum_right);

                for (size_t i = 0; i < AS; i++) {
                        carry[i] += (1.0-pm_right)*carry_left [i]*carry_right[AS];
                        carry[i] +=      pm_right *carry_left [i]*poly_sum_right;
                        carry[i] += (1.0-pm_right)*carry_right[i]*carry_left [AS];
                        carry[i] += (1.0-pm_right)*carry_left [i]*carry_right[i];
                }
                return carry;
        }
}

template <size_t AS, typename PC = double>
class pt_polynomial_t : public polynomial_t<AS, PC> {
public:
        pt_polynomial_t(const polynomial_t<AS, PC>& polynomial)
                : polynomial_t<AS, PC>(polynomial)
                { }
};

template <size_t AS, typename AC, typename PC>
pt_polynomial_t<AS, PC>
pt_likelihood(const pt_root_t& root, const std::vector<AC>& observations) {
        return tfbayes_detail::poly_sum<AS, AC, PC>(
                tfbayes_detail::likelihood_root<AS, AC, PC>(root, observations));
}

#endif /* __TFBAYES_PHYLOTREE_POLYNOMIAL_HH__ */
