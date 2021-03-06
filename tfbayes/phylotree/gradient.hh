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

#ifndef __TFBAYES_PHYLOTREE_GRADIENT_HH__
#define __TFBAYES_PHYLOTREE_GRADIENT_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm> /* std::min */
#include <cmath>
#include <vector>
#include <set>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/polynomial.hh>

namespace tfbayes_detail {
        template <size_t AS, typename PC>
        class derivatives_t : public std::vector<polynomial_t<AS, PC> > {
                derivatives_t()
                        : std::vector<polynomial_t<AS, PC> >()
                        { }
                derivatives_t(size_t k)
                        : std::vector<polynomial_t<AS, PC> >(k)
                        { }
        };
        template <size_t AS, typename AC, typename PC>
        class partial_t : public carry_t<AS, AC, PC> {
        public:
                partial_t(size_t n)
                        : derivatives(n)
                        { }

                std::vector<carry_t<AS, AC, PC> > derivatives;
        };
        template <size_t AS, typename AC, typename PC>
        partial_t<AS, AC, PC> derivative_leaf(
                const pt_leaf_t& leaf,
                const std::vector<AC>& observations,
                pt_node_t::id_t n) {
                partial_t<AS, AC, PC> partial(n);
                if (observations[leaf.id] == -1) {
                        // this leaf should be ignored, no nucleotide
                        // is present (missing data)
                        partial[AS] = 1.0;
                }
                else {
                        // leaf is present and has a nucleotide
                        partial[observations[leaf.id]] += 1.0;
                }
                return partial;
        }
        template <size_t AS, typename AC, typename PC>
        partial_t<AS, AC, PC> derivative_node(
                const pt_node_t& node,
                partial_t<AS, AC, PC>& partial_left,
                partial_t<AS, AC, PC>& partial_right,
                const std::vector<AC>& observations,
                pt_node_t::id_t n) {

                double  pm_left  = 1.0-std::exp(-node.left ().d);
                double  pm_right = 1.0-std::exp(-node.right().d);
                double  pn_left  = 1.0-pm_left;
                double  pn_right = 1.0-pm_right;

                double dpm_left  =  pn_left;
                double dpm_right =  pn_right;
                double dpn_left  = -pn_left;
                double dpn_right = -pn_right;

                partial_t<AS, AC, PC> partial(n);
                const polynomial_t<AS, PC> poly_sum_left  = poly_sum(partial_left);
                const polynomial_t<AS, PC> poly_sum_right = poly_sum(partial_right);

                partial[AS] +=
                        (pn_left *partial_left [AS] + pm_left *poly_sum_left)*
                        (pn_right*partial_right[AS] + pm_right*poly_sum_right);

                for (size_t i = 0; i < AS; i++) {
                        partial[i] += pn_left *pn_right*partial_left [i]*partial_right[AS];
                        partial[i] += pn_left *pm_right*partial_left [i]*poly_sum_right;
                        partial[i] += pn_right*pn_left *partial_right[i]*partial_left [AS];
                        partial[i] += pn_right*pm_left *partial_right[i]*poly_sum_left;
                        partial[i] += pn_left *pn_right*partial_left [i]*partial_right[i];
                }
                for (pt_node_t::id_t id = 0; id < n; id++)
                {
                        const polynomial_t<AS, PC> deri_sum_left  = poly_sum(partial_left .derivatives[id]);
                        const polynomial_t<AS, PC> deri_sum_right = poly_sum(partial_right.derivatives[id]);
                        /* Derivative of sigma
                         */
                        if (node.left().id == id) {
                                partial.derivatives[id][AS] +=
                                        (dpn_left*partial_left [AS]  + dpm_left *poly_sum_left)*
                                        (pn_right*partial_right[AS]  +  pm_right*poly_sum_right);
                        }
                        else if (node.right().id == id) {
                                partial.derivatives[id][AS] +=
                                        ( pn_left *partial_left [AS] +  pm_left *poly_sum_left)*
                                        (dpn_right*partial_right[AS] + dpm_right*poly_sum_right);
                        }
                        else {
                                partial.derivatives[id][AS] +=
                                        (pn_left *partial_left [AS]                 + pm_left *poly_sum_left)*
                                        (pn_right*partial_right.derivatives[id][AS] + pm_right*deri_sum_right);
                                partial.derivatives[id][AS] +=
                                        (pn_left *partial_left.derivatives [id][AS] + pm_left *deri_sum_left)*
                                        (pn_right*partial_right[AS]                 + pm_right*poly_sum_right);
                        }
                        /* Derivative of phi
                         */
                        if (node.left().id == id) {
                                for (size_t i = 0; i < AS; i++) {
                                        partial.derivatives[id][i] += dpn_left*pn_right*partial_left[AS]*partial_right[ i];
                                        partial.derivatives[id][i] += dpm_left*pn_right*poly_sum_left   *partial_right[ i];
                                        partial.derivatives[id][i] += dpn_left*pn_right*partial_left[ i]*partial_right[AS];
                                        partial.derivatives[id][i] += dpn_left*pm_right*partial_left[ i]*poly_sum_right;
                                        partial.derivatives[id][i] += dpn_left*pn_right*partial_left[ i]*partial_right[ i];
                                }
                        }
                        else if (node.right().id == id) {
                                for (size_t i = 0; i < AS; i++) {
                                        partial.derivatives[id][i] +=  pn_left*dpn_right*partial_left[ i]*partial_right[AS];
                                        partial.derivatives[id][i] +=  pn_left*dpm_right*partial_left[ i]*poly_sum_right;
                                        partial.derivatives[id][i] +=  pn_left*dpn_right*partial_left[AS]*partial_right[ i];
                                        partial.derivatives[id][i] +=  pm_left*dpn_right*poly_sum_left   *partial_right[ i];
                                        partial.derivatives[id][i] +=  pn_left*dpn_right*partial_left[ i]*partial_right[ i];
                                }
                        }
                        else {
                                for (size_t i = 0; i < AS; i++) {
                                        partial.derivatives[id][i] += pn_left *pn_right*partial_left[i]*partial_right.derivatives[id][AS];
                                        partial.derivatives[id][i] += pn_left *pm_right*partial_left[i]*deri_sum_right;

                                        partial.derivatives[id][i] += pn_right*pn_left*partial_right[i]*partial_left .derivatives[id][AS];
                                        partial.derivatives[id][i] += pn_right*pm_left*partial_right[i]*deri_sum_left;

                                        partial.derivatives[id][i] += pn_right*pn_left*partial_right.derivatives[id][i]*partial_left[AS];
                                        partial.derivatives[id][i] += pn_right*pm_left*partial_right.derivatives[id][i]*poly_sum_left;

                                        partial.derivatives[id][i] += pn_left *pn_right*partial_left.derivatives[id][i]*partial_right[AS];
                                        partial.derivatives[id][i] += pn_left *pm_right*partial_left.derivatives[id][i]*poly_sum_right;

                                        partial.derivatives[id][i] += pn_left *pn_right*partial_left [i]*partial_right.derivatives[id][i];
                                        partial.derivatives[id][i] += pn_left *pn_right*partial_right[i]*partial_left .derivatives[id][i];
                                }
                        }
                }
                return partial;
        }
        template <size_t AS, typename AC, typename PC>
        partial_t<AS, AC, PC> derivative_rec(
                const pt_node_t& node,
                const std::vector<AC>& observations,
                pt_node_t::id_t n) {
                if (node.leaf()) {
                        return derivative_leaf<AS, AC, PC>(static_cast<const pt_leaf_t&>(node), observations, n);
                }
                else {
                        partial_t<AS, AC, PC> partial_left  = derivative_rec<AS, AC, PC>(node.left (), observations, n);
                        partial_t<AS, AC, PC> partial_right = derivative_rec<AS, AC, PC>(node.right(), observations, n);

                        return derivative_node<AS, AC, PC>(node, partial_left, partial_right, observations, n);
                }
        }
        template <size_t AS, typename AC, typename PC>
        partial_t<AS, AC, PC> derivative_root(
                const pt_root_t& root,
                const std::vector<AC>& observations) {

                pt_node_t::id_t n = root.n_nodes-1;

                // if the root has no outgroup then treat it as a
                // simple node
                if (!root.outgroup()) {
                        derivative_rec<AS, AC, PC>(root, observations, n);
                }
                partial_t<AS, AC, PC> partial_left  = derivative_rec<AS, AC, PC>( root,            observations, n);
                partial_t<AS, AC, PC> partial_right = derivative_rec<AS, AC, PC>(*root.outgroup(), observations, n);

                double  pm_right = 1.0-std::exp(-root.outgroup()->d);
                double  pn_right = 1.0-pm_right;

                double dpm_right =  pn_right;
                double dpn_right = -pn_right;

                partial_t<AS, AC, PC> partial(n);
                const polynomial_t<AS, PC> poly_sum_left  = poly_sum(partial_left);
                const polynomial_t<AS, PC> poly_sum_right = poly_sum(partial_right);

                partial[AS] +=
                        partial_left [AS]*(pn_right*partial_right[AS] + pm_right*poly_sum_right);

                for (size_t i = 0; i < AS; i++) {
                        partial[i] += pn_right*partial_left [i]*partial_right[AS];
                        partial[i] += pm_right*partial_left [i]*poly_sum_right;
                        partial[i] += pn_right*partial_right[i]*partial_left [AS];
                        partial[i] += pn_right*partial_left [i]*partial_right[i];
                }
                for (pt_node_t::id_t id = 0; id < n; id++)
                {
                        const polynomial_t<AS, PC> deri_sum_left  = poly_sum(partial_left .derivatives[id]);
                        const polynomial_t<AS, PC> deri_sum_right = poly_sum(partial_right.derivatives[id]);
                        /* Derivative of sigma
                         */
                        if (root.outgroup()->id == id) {
                                partial.derivatives[id][AS] +=
                                                   partial_left [AS]*
                                        (dpn_right*partial_right[AS] + dpm_right*poly_sum_right);
                        }
                        else {
                                partial.derivatives[id][AS] +=
                                                  partial_left [AS]*
                                        (pn_right*partial_right.derivatives[id][AS] + pm_right*deri_sum_right);
                                partial.derivatives[id][AS] +=
                                                  partial_left.derivatives [id][AS]*
                                        (pn_right*partial_right[AS]                 + pm_right*poly_sum_right);
                        }
                        /* Derivative of phi
                         */
                        if (root.outgroup()->id == id) {
                                for (size_t i = 0; i < AS; i++) {
                                        partial.derivatives[id][i] +=  dpn_right*partial_left[ i]*partial_right[AS];
                                        partial.derivatives[id][i] +=  dpm_right*partial_left[ i]*poly_sum_right;
                                        partial.derivatives[id][i] +=  dpn_right*partial_left[AS]*partial_right[ i];
                                        partial.derivatives[id][i] +=  dpn_right*partial_left[ i]*partial_right[ i];
                                }
                        }
                        else {
                                for (size_t i = 0; i < AS; i++) {
                                        partial.derivatives[id][i] += pn_right*partial_left[i]*partial_right.derivatives[id][AS];
                                        partial.derivatives[id][i] += pm_right*partial_left[i]*deri_sum_right;

                                        partial.derivatives[id][i] += pn_right*partial_right[i]*partial_left .derivatives[id][AS];

                                        partial.derivatives[id][i] += pn_right*partial_right.derivatives[id][i]*partial_left[AS];

                                        partial.derivatives[id][i] += pn_right*partial_left.derivatives[id][i]*partial_right[AS];
                                        partial.derivatives[id][i] += pm_right*partial_left.derivatives[id][i]*poly_sum_right;

                                        partial.derivatives[id][i] += pn_right*partial_left [i]*partial_right.derivatives[id][i];
                                        partial.derivatives[id][i] += pn_right*partial_right[i]*partial_left .derivatives[id][i];
                                }
                        }
                }
                return partial;
        }
}

template <size_t AS, typename PC = double>
class pt_polynomial_derivative_t : public polynomial_t<AS, PC> {
public:
        typedef polynomial_t<AS, PC> base_t;

        pt_polynomial_derivative_t()
                : base_t      (),
                  _derivative_()
                { }
        pt_polynomial_derivative_t(size_t k)
                : base_t      (),
                  _derivative_(k)
                { }

        pt_polynomial_derivative_t operator=(const polynomial_t<AS, PC>& polynomial) {
                base_t::operator=(polynomial);
                return *this;
        }
        const std::vector<polynomial_t<AS, PC> >& derivative() const {
                return _derivative_;
        }
        std::vector<polynomial_t<AS, PC> >& derivative() {
                return _derivative_;
        }
        size_t d() const {
                return _derivative_.size();
        }
protected:
        std::vector<polynomial_t<AS, PC> > _derivative_;
};

class pt_marginal_derivative_t : std::pair<double, std::vector<double> > {
public:
        typedef std::pair<double, std::vector<double> > base_t;

        pt_marginal_derivative_t()
                : base_t()
                { }
        pt_marginal_derivative_t(size_t k)
                : base_t(0.0, std::vector<double>(k))
                { }
        pt_marginal_derivative_t(size_t k, double init)
                : base_t(init, std::vector<double>(k, init))
                { }

        pt_marginal_derivative_t& operator=(double value) {
                base_t::first = value;
                return *this;
        }
        operator double&() {
                return base_t::first;
        }
        const std::vector<double>& derivative() const {
                return base_t::second;
        }
        std::vector<double>& derivative() {
                return base_t::second;
        }
        size_t d() const {
                return base_t::second.size();
        }
};

template <size_t AS, typename AC, typename PC>
pt_polynomial_derivative_t<AS, PC>
pt_likelihood_derivative(const pt_root_t& root, const std::vector<AC>& observations)
{
        pt_node_t::id_t n = root.n_nodes-1;

        pt_polynomial_derivative_t<AS, PC> likelihood(n);
        
        tfbayes_detail::partial_t<AS, AC, PC> partial = tfbayes_detail::derivative_root<AS, AC, PC>(
                root, observations);

        likelihood = tfbayes_detail::poly_sum(partial);
 
        for (pt_node_t::id_t id = 0; id < n; id++) {
                likelihood.derivative()[id] = tfbayes_detail::poly_sum(partial.derivatives[id]);
        }
        return likelihood;
}

#endif /* __TFBAYES_PHYLOTREE_GRADIENT_HH__ */
