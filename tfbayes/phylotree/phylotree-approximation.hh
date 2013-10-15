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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_psi.h>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/phylotree/posterior.hh>

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
struct kl_divergence_data {
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& variational;
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& likelihood;
        const exponent_t<ALPHABET_SIZE, CODE_TYPE>& alpha;
        boost::array<double, ALPHABET_SIZE>& theta;
};

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
double
dkl_f(double * x, size_t dim, void * params)
{
        struct kl_divergence_data<ALPHABET_SIZE, CODE_TYPE>* data =
                (kl_divergence_data<ALPHABET_SIZE, CODE_TYPE> *)(params);

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
                double p = pt_posterior_density<ALPHABET_SIZE, CODE_TYPE>(
                        data->likelihood, data->alpha, data->theta);
                double q = pt_posterior_density<ALPHABET_SIZE, CODE_TYPE>(
                        data->variational, data->alpha, data->theta);

                return q*log(p);
        }
}

/* Compute the entropy of the variational distribution. */
template <size_t ALPHABET_SIZE, typename CODE_TYPE>
double
dkl_variational_entropy(
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& variational,
        const exponent_t<ALPHABET_SIZE, CODE_TYPE>& alpha)
{
        assert(variational.size() == 1);

        /* fetch the exponent of the variational distribution */
        const exponent_t<ALPHABET_SIZE, CODE_TYPE>& exponent = variational.begin()->exponent();

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

/* Compute the Kullback-Leibler divergence * D(q||p), where q is the
 * variational distribution and p the actual distribution. This
 * function uses numerical integration to solve the integral. The
 * entropy of q can be solved analytically.
 */
template <size_t ALPHABET_SIZE, typename CODE_TYPE>
double
dkl(
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& variational,
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& likelihood,
        const exponent_t<ALPHABET_SIZE, CODE_TYPE>& alpha)
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
                dkl_f<ALPHABET_SIZE, CODE_TYPE>;

        struct kl_divergence_data<ALPHABET_SIZE, CODE_TYPE> data = {
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
        result += dkl_variational_entropy<ALPHABET_SIZE, CODE_TYPE>(variational, alpha);

        return -result;
}

/* Define a line that goes through the origin and the approximated
 * minimum (lower bound) of the Kullback-Leibler divergence. The
 * actual minimum is assumed to be somewhere on this line close to
 * the approximated point.
 */
template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE>
dkl_line(
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& variational,
        const exponent_t<ALPHABET_SIZE, CODE_TYPE>& alpha,
        const double lambda)
{
        polynomial_term_t<ALPHABET_SIZE, CODE_TYPE> term(*variational.begin());

        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                term.exponent()[i] += alpha[i];
                term.exponent()[i] *= lambda;
                term.exponent()[i] -= alpha[i];
        }
        return term;
}

/* Perform a simple line search to find the actual minimum of the
 * Kullback-Leibler divergence.
 */
template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE>
dkl_line_search(
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& variational,
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& likelihood,
        const exponent_t<ALPHABET_SIZE, CODE_TYPE>& alpha,
        const size_t n)
{
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> result = dkl_line<ALPHABET_SIZE, CODE_TYPE>(variational, alpha, 1.0);
        /* step size */
        double epsilon = 0.01;
        /* how much to scale the step size */
        double eta1    = 1.05;
        double eta2    = 0.50;
        /* position on the line */
        double lambda  = 1.0;
        double kl = dkl<ALPHABET_SIZE, CODE_TYPE>(result, likelihood, alpha);
        double kl_new;

        for (size_t i = 0; i < n; i++) {
                result = dkl_line<ALPHABET_SIZE, CODE_TYPE>(variational, alpha, lambda);
                kl_new = dkl<ALPHABET_SIZE, CODE_TYPE>(result, likelihood, alpha);
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

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE>
dkl_approximate(
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& poly)
{
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> result;
        polynomial_term_t<ALPHABET_SIZE, CODE_TYPE> result_term(1.0);
        double norm = 0;

        /* compute normalization constant */
        for (typename polynomial_t<ALPHABET_SIZE, CODE_TYPE>::const_iterator it = poly.begin(); it != poly.end(); it++)
        {
                norm += it->coefficient();
        }
        /* compute approximation */
        for (typename polynomial_t<ALPHABET_SIZE, CODE_TYPE>::const_iterator it = poly.begin(); it != poly.end(); it++)
        {
                /* loop over the alphabet */
                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        result_term.exponent()[i] += it->coefficient()/norm * it->exponent()[i];
                }
        }
        result += result_term;

        return result;
}

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
boost::array<double, ALPHABET_SIZE>
dkl_optimize_params(
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& poly)
{
        boost::array<double, ALPHABET_SIZE> result;
        double norm = -HUGE_VAL;

        /* initialize result */
        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                result[i] = -HUGE_VAL;
        }

        /* compute normalization constant */
        for (typename polynomial_t<ALPHABET_SIZE, CODE_TYPE>::const_iterator it = poly.begin(); it != poly.end(); it++)
        {
                norm = logadd(norm, it->coefficient() + mbeta_log(it->exponent()));
        }
        /* for each member of the alphabet */
        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                /* compute approximation */
                for (typename polynomial_t<ALPHABET_SIZE, CODE_TYPE>::const_iterator it = poly.begin(); it != poly.end(); it++)
                {
                        /* pi Beta(alpha_i) */
                        double tmp1 = it->coefficient() + mbeta_log(it->exponent());
                        double sum  = 0;
                        /* sum the exponents */
                        for (size_t j = 0; j < ALPHABET_SIZE; j++) {
                                sum += it->exponent()[j];
                        }
                        /* digamma(alpha_i,k) - digamma(alpha_i,0),
                         * but since this is negative we need to move
                         * the negation outside the sum */
                        double tmp2 = log(gsl_sf_psi(sum) - gsl_sf_psi(it->exponent()[i]));
                        /* multiply the two partial results */
                        result[i] = logadd(result[i], tmp1 + tmp2);
                }
        }
        /* normalize the result and leave log scale */
        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                /* negate here */
                result[i] = -exp(result[i] - norm);
        }
        return result;
}


template <size_t ALPHABET_SIZE, typename CODE_TYPE>
int
dkl_optimize_f(
        const gsl_vector * xi,
        void *params,
        gsl_vector * f)
{
        boost::array<double, ALPHABET_SIZE>& p =
                *(boost::array<double, ALPHABET_SIZE>*)params;

        double sum = 0;

        /* check position */
        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                if (gsl_vector_get(xi, i) < 0.0) {
                        return GSL_FAILURE;
                }
        }

        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                sum += gsl_vector_get(xi, i);
        }

        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                /* psi(xi_i) - psi(xi_0) - p_i */
                const double xi_i = gsl_vector_get(xi, i);
                const double  f_i =
                        gsl_sf_psi(xi_i) - gsl_sf_psi(sum) - p[i];
                gsl_vector_set(f, i, f_i);
        }
     
        return GSL_SUCCESS;
}

/* Optimize the variational distribution q such that the
 * Kullback-Leibler divergence D(p||q) is minimized. This function
 * uses a simple root-finding algorithm to find the optimum.
 */
template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE>
dkl_optimize(
        const polynomial_t<ALPHABET_SIZE, CODE_TYPE>& _poly,
        const exponent_t<ALPHABET_SIZE, CODE_TYPE>& alpha)
{
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> poly = _poly.normalize();

        /* integrate alpha into the polynomial */
        poly *= alpha;

        /* use the approximation as an initial value for optimizing
         * the Kullback-Leibler divergence */
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> variational =
                dkl_approximate<ALPHABET_SIZE, CODE_TYPE>(poly);

        /* assert that this really is a variational distribution with
         * only one component */
        assert(variational.size() == 1);

        /* copy the initial value to a GSL vector */
        gsl_vector *x = gsl_vector_alloc(ALPHABET_SIZE);

        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                gsl_vector_set(x, i, variational.begin()->exponent()[i]);
        }

        /* the parameters of the objective function */
        boost::array<double, ALPHABET_SIZE> params =
                dkl_optimize_params(poly);

        /* the root finding algorithm from GSL */
        const gsl_multiroot_fsolver_type *T;
        gsl_multiroot_fsolver *s;

        /* some auxiliary variables */
        int status;
        size_t iter = 0;
     
        gsl_multiroot_function f = {
                /* the objective function */
                dkl_optimize_f<ALPHABET_SIZE, CODE_TYPE>,
                /* dimensionality of the problem */
                ALPHABET_SIZE,
                /* parameters for the objective function */
                &params
        };
          
        T = gsl_multiroot_fsolver_hybrids;
        s = gsl_multiroot_fsolver_alloc(T, ALPHABET_SIZE);
        gsl_multiroot_fsolver_set(s, &f, x);

        do {
                iter++;
                status = gsl_multiroot_fsolver_iterate (s);

                if (status)   /* check if solver is stuck */
                        break;

                status = gsl_multiroot_test_residual(s->f, 1e-7);
        }
        while (status == GSL_CONTINUE && iter < 1000);

        /* copy the result to a polynomial */
        exponent_t<ALPHABET_SIZE, CODE_TYPE> result;

        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                result[i] = gsl_vector_get(x, i);
        }

        /* free GSL data structures */
        gsl_multiroot_fsolver_free(s);
        gsl_vector_free(x);

        return result;
}

#endif /* PHYLOTREE_APPROXIMATION_HH */
