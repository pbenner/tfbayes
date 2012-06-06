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

#ifndef PHYLOTREE_SAMPLER_HH
#define PHYLOTREE_SAMPLER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cmath>
#include <algorithm> /* std::min */
#include <iomanip>

#include <tfbayes/exception.h>
#include <tfbayes/polynomial.hh>

#include <alignment.hh>
#include <phylotree.hh>
#include <phylotree-distribution.hh>
#include <phylotree-gradient-coefficient.hh>
#include <clonable.hh>

#include <sys/time.h>

#include <gsl/gsl_randist.h>

class jumping_distribution_t : public clonable
{
public:
        jumping_distribution_t* clone() const = 0;

        virtual double p(double d_old, double d_new) const = 0;
        virtual double sample(gsl_rng * rng, double d_old) const = 0;
        virtual void increase_jump(double eta) = 0;
        virtual void decrease_jump(double eta) = 0;
};

class normal_jump_t : public jumping_distribution_t
{
public:
        normal_jump_t(double sigma_square = 0.05)
                : sigma_square(sigma_square) { }

        normal_jump_t* clone() const {
                return new normal_jump_t(*this);
        }

        double p(double d_old, double d_new) const {
                return 1.0;
        }
        double sample(gsl_rng * rng, double d_old) const {
                return d_old+gsl_ran_gaussian(rng, sigma_square);
        }
        void increase_jump(double eta) {
                sigma_square = sigma_square+eta;
        }
        void decrease_jump(double eta) {
                sigma_square = std::max(sigma_square-eta, 0.0);
        }
private:
        double sigma_square;
};

class gamma_jump_t : public jumping_distribution_t
{
public:
        gamma_jump_t(double r, double lambda)
                : gamma_distribution(r, lambda) { }

        gamma_jump_t* clone() const {
                return new gamma_jump_t(*this);
        }

        // p(d_new -> d_old)/p(d_old -> d_new) = p(d_old)/p(d_new)
        double p(double d_old, double d_new) const {
                
                return gamma_distribution.pdf(d_old)/gamma_distribution.pdf(d_new);
        }
        double sample(gsl_rng * rng, double d_old) const {
                return gamma_distribution.sample(rng);
        }
        void increase_jump(double eta) {
        }
        void decrease_jump(double eta) {
        }
private:
        gamma_distribution_t gamma_distribution;
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_metropolis_hastings_t
{
public:
        pt_metropolis_hastings_t(alignment_t<CODE_TYPE>& alignment,
                                 const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
                                 double r, double lambda,
                                 jumping_distribution_t& jumping_distribution,
                                 double acceptance_rate = 0.7)
                : alignment(alignment), alpha(alpha),
                  gamma_distribution(r, lambda),
                  acceptance_rate(acceptance_rate) {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize nodes and jumping distributions
                nodes = alignment.tree->get_nodes();
                for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                        jumping_distributions[*is] = jumping_distribution.clone();
                }
        }

        ~pt_metropolis_hastings_t() {
                // free random generator
                gsl_rng_free(rng);
                // free jumping distributions
                for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                        delete(jumping_distributions[*is]);
                }
        }

        double log_likelihood() {
                double result = 0;
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        const pt_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(alignment.tree);
                        // loop over monomials
                        double tmp = -HUGE_VAL;
                        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator ut = polynomial.begin();
                             ut != polynomial.end(); ut++) {
                                tmp = logadd(tmp, log(ut->coefficient()) + mbeta_log(ut->exponent(), alpha) - mbeta_log(alpha));
                        }
                        result += tmp;
                }
                return result;
        }
        std::pair<double,bool> sample(pt_node_t* node, double log_likelihood_ref) {
                double rho;
                double x;
                double log_likelihood_new;

                double d_old = node->d;
                double d_new = jumping_distributions[node]->sample(rng, d_old);

                // compute new log likelihood
                node->d            = d_new;
                log_likelihood_new = log_likelihood();

                // compute acceptance probability
                rho = exp(log_likelihood_new-log_likelihood_ref)
                        *gamma_distribution.pdf(d_new)/gamma_distribution.pdf(d_old)
                        *jumping_distributions[node]->p(d_old, d_new);
                x   = gsl_ran_flat(rng, 0.0, 1.0);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        print_debug("accepted: %f\n", d_new);
                        return std::pair<double, bool>(log_likelihood_new, true);
                }
                else {
                        // sample rejected
                        print_debug("rejected: %f\n", d_new);
                        node->d = d_old;
                        return std::pair<double, bool>(log_likelihood_ref, false);
                }
        }
        double sample(boost::unordered_map<pt_node_t*, double>& means,
                      boost::unordered_map<pt_node_t*, double>& acceptance,
                      size_t i, bool burnin) {
                double log_likelihood_ref = log_likelihood();
                // loop over nodes
                for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                        std::pair<double, bool> result = sample(*is, log_likelihood_ref);
                        log_likelihood_ref = result.first;
                        // estimate acceptance rate
                        if (result.second) {
                                acceptance[*is] = (acceptance[*is]*i + 1.0)/(i+1.0);
                        }
                        else {
                                acceptance[*is] = (acceptance[*is]*i)/(i+1.0);
                        }
                }
                // update means
                for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                        means[*is] = ((double)i*means[*is] + (*is)->d)/(double)(i+1);
                        print_debug("rate: %f\n", acceptance[*is]);
                        if (burnin) {
                                if (acceptance[*is] > acceptance_rate) {
                                        jumping_distributions[*is]->increase_jump(fabs(acceptance[*is]-acceptance_rate));
                                }
                                else {
                                        jumping_distributions[*is]->decrease_jump(fabs(acceptance[*is]-acceptance_rate));
                                }
                        }
                        std::cerr << std::setprecision(8)
                                  << std::fixed
                                  << means[*is] << " ";
                }
                std::cerr << std::endl;

                return log_likelihood_ref;
        }
        double sample(size_t burnin, size_t n) {
                boost::unordered_map<pt_node_t*, double> means;
                boost::unordered_map<pt_node_t*, double> acceptance;
                double result;
                // burn in
                for (size_t i = 0; i < burnin; i++) {
                        sample(means, acceptance, i, true);
                }
                // sample n times
                for (size_t i = 0; i < n; i++) {
                        result = sample(means, acceptance, i, false);
                }
                // insert means into the tree
                for (boost::unordered_map<pt_node_t*, double>::iterator it = means.begin(); it != means.end(); it++) {
                        it->first->d = it->second;
                }
                return result;
        }

private:
        alignment_t<CODE_TYPE>& alignment;
        exponent_t<CODE_TYPE, ALPHABET_SIZE> alpha;

        gamma_distribution_t gamma_distribution;

        pt_node_t::nodes_t nodes;
        boost::unordered_map<pt_node_t*, jumping_distribution_t*> jumping_distributions;
        double acceptance_rate;

        gsl_rng * rng;
};

#endif /* PHYLOTREE_SAMPLER_HH */
