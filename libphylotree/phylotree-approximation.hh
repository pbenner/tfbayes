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
#include <posterior.hh>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_psi.h>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
struct kl_divergence_data {
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& variational;
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& likelihood;
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha;
        boost::array<double, ALPHABET_SIZE>& theta;
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
double kl_divergence_f(double * x, size_t dim, void * params)
{
        struct kl_divergence_data<CODE_TYPE, ALPHABET_SIZE>* data =
                (kl_divergence_data<CODE_TYPE, ALPHABET_SIZE> *)(params);

        data->theta[ALPHABET_SIZE-1] = 1.0;
        for (size_t i = 0; i < ALPHABET_SIZE-1; i++) {
                data->theta[i]                = x[i];
                data->theta[ALPHABET_SIZE-1] -= x[i];
        }

        /* check if we are inside the simplex */
        if (data->theta[ALPHABET_SIZE-1] < 0.0) {
                return 0.0;
        }
        else {
                double p = pt_posterior_density<CODE_TYPE, ALPHABET_SIZE>(
                        data->likelihood, data->alpha, data->theta);
                double q = pt_posterior_density<CODE_TYPE, ALPHABET_SIZE>(
                        data->variational, data->alpha, data->theta);

                return q*log(p);
        }
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
double variational_entropy(
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& variational,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha)
{
        /* fetch the exponent of the variational distribution */
        assert(variational.size() == 1);
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& exponent = variational.begin()->exponent();

        double sum    = 0.0;
        double result = mbeta_log(exponent, alpha);

        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                double tmp = exponent[i];
                result -= (tmp + alpha[i] - 1.0)*gsl_sf_psi(tmp + alpha[i]);
                sum    +=  tmp + alpha[i];
        }
        result += (sum - ALPHABET_SIZE)*gsl_sf_psi(sum);

        return result;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
double kl_divergence(
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& variational,
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& likelihood,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha)
{
        boost::array<double, ALPHABET_SIZE> theta;
        double xl[ALPHABET_SIZE-1];
        double xu[ALPHABET_SIZE-1];
        const gsl_rng_type *T;
        gsl_rng *r;

        /* compute Int Q Log P dw                                              *
         ***********************************************************************/
        size_t calls = 500000;
        double result, err;

        gsl_monte_function F;

        /* generate the function that we want to integrate */
        double (*f)(double *, size_t, void *) =
                kl_divergence_f<CODE_TYPE, ALPHABET_SIZE>;

        struct kl_divergence_data<CODE_TYPE, ALPHABET_SIZE> data = {
                variational, likelihood, alpha, theta
        };

        for (size_t i = 0; i < ALPHABET_SIZE-1; i++) {
                xl[i] = 0.0;
                xu[i] = 1.0;
        }

        F.f      = f;
        F.dim    = ALPHABET_SIZE-1;
        F.params = &data;
     
        gsl_rng_env_setup ();
     
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);

        gsl_monte_miser_state *s = gsl_monte_miser_alloc (ALPHABET_SIZE-1);
        gsl_monte_miser_integrate (&F, xl, xu, ALPHABET_SIZE-1, calls, r, s, 
                                    &result, &err);
        gsl_monte_miser_free (s);

        /* add - Int Q Log Q dw                                                *
         ***********************************************************************/
        result += variational_entropy<CODE_TYPE, ALPHABET_SIZE>(variational, alpha);

        return -result;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_line(
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& variational,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
        const double lambda)
{
        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> term(*variational.begin());

        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                term.exponent()[i] += alpha[i];
                term.exponent()[i] *= lambda;
                term.exponent()[i] -= alpha[i];
        }
        return term;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_line_search(
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& variational,
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& likelihood,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
        const size_t n)
{
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> result = pt_line<CODE_TYPE, ALPHABET_SIZE>(variational, alpha, 1.0);
        /* step size */
        double epsilon = 0.01;
        /* how much to scale the step size */
        double eta1    = 1.05;
        double eta2    = 0.50;
        /* position on the line */
        double lambda  = 1.0;
        double kl = kl_divergence<CODE_TYPE, ALPHABET_SIZE>(result, likelihood, alpha);
        double kl_new;

        for (size_t i = 0; i < n; i++) {
                result = pt_line<CODE_TYPE, ALPHABET_SIZE>(variational, alpha, lambda);
                kl_new = kl_divergence<CODE_TYPE, ALPHABET_SIZE>(result, likelihood, alpha);
                if (kl_new < kl) {
                        /* go faster */
                        epsilon *= eta1;
                }
                else {
                        /* go slower */
                        epsilon *= eta2;
                        /* turn around */
                        epsilon *= -1;
                }
                lambda = (1.0-epsilon)*lambda;
                kl     = kl_new;
        }

        return result;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_approximate(
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& poly)
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
