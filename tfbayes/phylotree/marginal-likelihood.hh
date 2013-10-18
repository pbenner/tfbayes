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

#ifndef MARGINAL_LIKELIHOOD_HH
#define MARGINAL_LIKELIHOOD_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <limits>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/utility/statistics.hh>
#include <tfbayes/utility/logarithmetic.hh>
#include <tfbayes/utility/polynomial.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

template <size_t AS, typename PC>
double pt_marginal_likelihood(
        const polynomial_t<AS, PC>& polynomial,
        const exponent_t<AS, PC>& alpha)
{
        double result = std::numeric_limits<double>::min();
        double mbeta_alpha = mbeta_log(alpha);

        for (typename polynomial_t<AS, PC>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                result = logadd(result, log(it->coefficient()) +
                                mbeta_log(it->exponent(), alpha) - mbeta_alpha);
        }

        return result;
}

template <size_t AS, typename AC, typename PC>
double pt_marginal_likelihood(
        const pt_root_t& node,
        const std::vector<AC>& observations,
        const exponent_t<AS, PC>& alpha)
{
        const polynomial_t<AS, PC> polynomial =
                pt_polynomial_t<AS, AC, PC>(node, observations);

        return pt_marginal_likelihood<AS, PC>(polynomial, alpha);
}

#endif /* MARGINAL_LIKELIHOOD_HH */
