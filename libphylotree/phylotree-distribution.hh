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

#ifndef PHYLOTREE_DISTRIBUTION_HH
#define PHYLOTREE_DISTRIBUTION_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>

class gamma_distribution_t {
public:
        gamma_distribution_t(double r, double lambda)
                : r(r), r_gamma(gsl_sf_gamma(r)), lambda(lambda),
                  lambda_pow_r(pow(lambda, r)) { }

        double pdf(double d) const {
                if (d > 0.0) {
                        return 1.0/lambda_pow_r * 1.0/r_gamma * pow(d, r-1.0) * exp(-d/lambda);
                }
                else {
                        return 0.0;
                }
        }
        double sample(const gsl_rng* rng) const {
                return gsl_ran_gamma(rng, r, lambda);
        }
        double log_gradient(double d) const {
                return (r-1.0)/d - 1.0/lambda;
        }

        double r;
        double r_gamma;
        double lambda;
        double lambda_pow_r;
};

#endif /* PHYLOTREE_DISTRIBUTION_HH */
