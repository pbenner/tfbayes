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
        virtual double sample(gsl_rng * rng, double d_old, double sigma) const = 0;
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
        double sample(gsl_rng * rng, double d_old, double sigma) const {
                return d_old+gsl_ran_gaussian(rng, sigma);
        }
        void increase_jump(double eta) {
                sigma_square = sigma_square+eta;
                print_debug("increasing sigma_square: %f\n", sigma_square);
        }
        void decrease_jump(double eta) {
                sigma_square = std::max(sigma_square-eta, 0.0);
                print_debug("decreasing sigma_square: %f\n", sigma_square);
        }
private:
        double sigma_square;
};

class geometric_jump_t : public jumping_distribution_t
{
public:
        geometric_jump_t(double sigma_square = 0.05)
                : sigma_square(sigma_square) { }

        geometric_jump_t* clone() const {
                return new geometric_jump_t(*this);
        }

        double p(double d_old, double d_new) const {
                return 1.0;
        }
        double sample(gsl_rng * rng, double d_old) const {
                return d_old+gsl_ran_gaussian(rng, sigma_square);
        }
        double sample(gsl_rng * rng, double d_old, double sigma) const {
                return d_old+gsl_ran_gaussian(rng, sigma);
        }
        void increase_jump(double eta) {
        }
        void decrease_jump(double eta) {
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
        double sample(gsl_rng * rng, double d_old, double sigma) const {
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
class pt_metropolis_hastings_t : public clonable
{
public:
        pt_metropolis_hastings_t(const pt_root_t* tree,
                                 const alignment_t<CODE_TYPE>& alignment,
                                 const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
                                 double r, double lambda,
                                 const jumping_distribution_t& jumping_distribution,
                                 double acceptance_rate = 0.7)
                : means(tree->size(),      0.0),
                  acceptance(tree->size(), 0.0),
                  history(tree->size(),    std::vector<double>()),
                  alignment(alignment),
                  alpha(alpha),
                  gamma_distribution(r, lambda),
                  acceptance_rate(acceptance_rate),
                  step(0),
                  tree(tree->clone()) {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize jumping distributions
                for (pt_node_t::id_t id = 0; id < tree->size(); id++) {
                        jumping_distributions.push_back(jumping_distribution.clone());
                }
        }
        pt_metropolis_hastings_t(const pt_metropolis_hastings_t& mh)
                : means(mh.means),
                  acceptance(mh.acceptance),
                  history(mh.history),
                  alignment(mh.alignment),
                  alpha(mh.alpha),
                  gamma_distribution(mh.gamma_distribution),
                  jumping_distributions(mh.jumping_distributions),
                  acceptance_rate(mh.acceptance_rate),
                  step(mh.step),
                  tree(mh.tree->clone())
                {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize nodes and jumping distributions
                for (pt_node_t::id_t id = 0; id < tree->size(); id++) {
                        jumping_distributions[id] = mh.jumping_distributions[id]->clone();
                }
        }
        ~pt_metropolis_hastings_t() {
                // free random generator
                gsl_rng_free(rng);
                // free jumping distributions
                for (pt_node_t::id_t id = 0; id < tree->size(); id++) {
                        delete(jumping_distributions[id]);
                }
                tree->destroy();
        }

        pt_metropolis_hastings_t* clone() const {
                return new pt_metropolis_hastings_t(*this);
        }

        double log_likelihood() {
                double result = 0;
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        it.apply(tree);
                        const pt_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(tree);
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
        void update_means(bool print) {
                // update root
                history[0].push_back(0.0);
                // update means
                for (pt_node_t::id_t id = 1; id < tree->size(); id++) {
                        pt_node_t* node = tree->node_map[id];
                        means[id] = ((double)step*means[id] + node->d)/(double)(step+1);
                        history[id].push_back(means[id]);
                        if (print) {
                                std::cout << std::setprecision(8)
                                          << std::fixed
                                          << means[id] << " ";
                        }
                }
                if (print) {
                        std::cout << std::endl;
                }
        }
        void update_steps() {
                for (pt_node_t::id_t id = 1; id < tree->size(); id++) {
                        if (acceptance[id] > acceptance_rate) {
                                jumping_distributions[id]->increase_jump(fabs(acceptance[id]-acceptance_rate));
                        }
                        else {
                                jumping_distributions[id]->decrease_jump(fabs(acceptance[id]-acceptance_rate));
                        }
                }
        }
        std::pair<double,bool> sample(pt_node_t::id_t id, double log_likelihood_ref) {
                double rho;
                double x;
                double log_likelihood_new;

                // generate a proposal
                double d_old = tree->node_map[id]->d;
                double d_new = std::min(0.0, jumping_distributions[id]->sample(rng, d_old));
                tree->node_map[id]->d = d_new;
                print_debug("old sample: %f\n", d_old);

                log_likelihood_new = log_likelihood();

                // compute acceptance probability
                rho = exp(log_likelihood_new-log_likelihood_ref)
                        *gamma_distribution.pdf(d_new)/gamma_distribution.pdf(d_old)
                        *jumping_distributions[id]->p(d_old, d_new);
                x   = gsl_ran_flat(rng, 0.0, 1.0);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        print_debug("accepted: %f\n", d_new);
                        return std::pair<double, bool>(log_likelihood_new, true);
                }
                else {
                        // sample rejected
                        tree->node_map[id]->d = d_old;
                        print_debug("rejected: %f\n", d_new);
                        return std::pair<double, bool>(log_likelihood_ref, false);
                }
        }
        void sample(bool print) {
                double log_likelihood_ref = log_likelihood();
                flockfile(stderr);
                std::cerr << "step: "
                          << step
                          << std::endl;
                fflush(stderr);
                funlockfile(stderr);
                // loop over nodes
                for (pt_node_t::id_t id = 1; id < tree->size(); id++) {
                        std::pair<double, bool> result = sample(id, log_likelihood_ref);
                        log_likelihood_ref = result.first;
                        // estimate acceptance rate
                        if (result.second) {
                                acceptance[id] = (acceptance[id]*step + 1.0)/(step+1.0);
                        }
                        else {
                                acceptance[id] = (acceptance[id]*step)/(step+1.0);
                        }
                }
                update_means(print);
                step++;
        }
        void burnin(size_t n, bool print = true) {
                // burn in
                for (size_t i = 0; i < n; i++) {
                        sample(print);
                        update_steps();
                }
                step = 0;
        }
        void sample(size_t n, bool print = true) {
                // sample n times
                for (size_t i = 0; i < n; i++) {
                        sample(print);
                }
        }
        void apply(pt_root_t* tree) {
                // insert means into the tree
                for (pt_root_t::node_map_t::const_iterator is = ++tree->node_map.begin(); is != tree->node_map.end(); is++) {
                        (*is)->d = means[(*is)->id];
                }
        }

        std::vector<double> means;
        std::vector<double> acceptance;
        std::vector<std::vector<double> > history;
private:
        const alignment_t<CODE_TYPE>& alignment;
        exponent_t<CODE_TYPE, ALPHABET_SIZE> alpha;

        gamma_distribution_t gamma_distribution;

        std::vector<jumping_distribution_t*> jumping_distributions;
        double acceptance_rate;

        gsl_rng * rng;
        size_t step;

        pt_root_t* tree;
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_geometric_hastings_t : public clonable
{
public:
        pt_geometric_hastings_t(const pt_root_t* tree,
                                const alignment_t<CODE_TYPE>& alignment,
                                const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
                                double r, double lambda,
                                const jumping_distribution_t& jumping_distribution,
                                double acceptance_rate = 0.7)
                : means(tree->size(),      0.0),
                  acceptance(tree->size(), 0.0),
                  history(tree->size(),    std::vector<double>()),
                  alignment(alignment),
                  alpha(alpha),
                  gamma_distribution(r, lambda),
                  acceptance_rate(acceptance_rate),
                  step(0),
                  tree(tree->clone()) {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize jumping distributions
                for (pt_node_t::id_t id = 0; id < tree->size(); id++) {
                        jumping_distributions.push_back(jumping_distribution.clone());
                }
        }
        pt_geometric_hastings_t(const pt_geometric_hastings_t& mh)
                : means(mh.means),
                  acceptance(mh.acceptance),
                  history(mh.history),
                  alignment(mh.alignment),
                  alpha(mh.alpha),
                  gamma_distribution(mh.gamma_distribution),
                  jumping_distributions(mh.jumping_distributions),
                  acceptance_rate(mh.acceptance_rate),
                  step(mh.step),
                  tree(mh.tree->clone())
                {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize nodes and jumping distributions
                for (pt_node_t::id_t id = 0; id < tree->size(); id++) {
                        jumping_distributions[id] = mh.jumping_distributions[id]->clone();
                }
        }
        ~pt_geometric_hastings_t() {
                // free random generator
                gsl_rng_free(rng);
                // free jumping distributions
                for (pt_node_t::id_t id = 0; id < tree->size(); id++) {
                        delete(jumping_distributions[id]);
                }
                // free tree
                tree->destroy();
        }

        pt_geometric_hastings_t* clone() const {
                return new pt_geometric_hastings_t(*this);
        }

        double log_likelihood() {
                double result = 0;
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        it.apply(tree);
                        const pt_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(tree);
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
        void update_means(bool print) {
                // update root
                history[0].push_back(0.0);
                // update means
                for (pt_node_t::id_t id = 1; id < tree->size(); id++) {
                        pt_node_t* node = tree->node_map[id];
                        means[id] = ((double)step*means[id] + node->d)/(double)(step+1);
                        history[id].push_back(means[id]);
                        if (print) {
                                std::cout << std::setprecision(8)
                                          << std::fixed
                                          << node->d << " ";
//                                          << means[id] << " ";
                        }
                }
                if (print) {
                        std::cout << std::endl;
                }
        }
        void update_steps() {
                for (pt_node_t::id_t id = 1; id < tree->size(); id++) {
                        if (acceptance[id] > acceptance_rate) {
                                jumping_distributions[id]->increase_jump(fabs(acceptance[id]-acceptance_rate));
                        }
                        else {
                                jumping_distributions[id]->decrease_jump(fabs(acceptance[id]-acceptance_rate));
                        }
                }
        }
        void update_acceptance(pt_node_t::id_t id, bool accepted) {
                // estimate acceptance rate
                if (accepted) {
                        acceptance[id] = (acceptance[id]*step + 1.0)/(step+1.0);
                }
                else {
                        acceptance[id] = (acceptance[id]*step)/(step+1.0);
                }
                print_debug("acceptance: %d=%f\n",  (int)id, acceptance[id]);
        }
        double sample_eta(pt_node_t::id_t id, double log_likelihood_ref) {
                double rho;
                double x;
                double log_likelihood_new;

                pt_node_t::id_t id_left  = tree->node_map[id]->left ->id;
                pt_node_t::id_t id_right = tree->node_map[id]->right->id;
                jumping_distribution_t* jumping_distribution =
                        jumping_distributions[id_left];

                // save old branch lengths
                double d1_old = tree->node_map[id_left ]->d;
                double d2_old = tree->node_map[id_right]->d;
                // generate a proposal
                double eta_old = (d1_old + d2_old)/2.0;
                double  xi_old = (d1_old - d2_old)/2.0;
                double eta_new = std::max(fabs(xi_old), jumping_distribution->sample(rng, eta_old, 0.01+eta_old/2.0));
                double  d1_new = eta_new + xi_old;
                double  d2_new = eta_new - xi_old;
                print_debug(" d1_old: %f\n",  d1_old);
                print_debug(" d2_old: %f\n",  d2_old);
                print_debug(" d1_new: %f\n",  d1_new);
                print_debug(" d2_new: %f\n",  d2_new);
                print_debug("eta_old: %f\n", eta_old);
                print_debug(" xi_old: %f\n",  xi_old);
                print_debug("eta_new: %f\n", eta_new);
                tree->node_map[id_left ]->d = d1_new;
                tree->node_map[id_right]->d = d2_new;

                log_likelihood_new = log_likelihood();

                rho = exp(log_likelihood_new-log_likelihood_ref)
                        *gamma_distribution.pdf(d1_new)/gamma_distribution.pdf(d1_old)
                        *gamma_distribution.pdf(d2_new)/gamma_distribution.pdf(d2_old)
                        *jumping_distribution->p(eta_old, eta_new);
                print_debug("    rho: %f\n",  rho);
                x = gsl_ran_flat(rng, 0.0, 1.0);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        print_debug("accepted d1 (eta): %f\n", d1_new);
                        print_debug("accepted d2 (eta): %f\n", d2_new);
                        update_acceptance(id_left, true);
                        return log_likelihood_new;
                }
                else {
                        // sample rejected
                        tree->node_map[id_left ]->d = d1_old;
                        tree->node_map[id_right]->d = d2_old;
                        print_debug("rejected d1 (eta): %f\n", d1_new);
                        print_debug("rejected d2 (eta): %f\n", d2_new);
                        update_acceptance(id_left, false);
                        return log_likelihood_ref;
                }
        }
        double sample_xi(pt_node_t::id_t id, double log_likelihood_ref) {
                double rho;
                double x;
                double log_likelihood_new;

                pt_node_t::id_t id_left  = tree->node_map[id]->left ->id;
                pt_node_t::id_t id_right = tree->node_map[id]->right->id;
                jumping_distribution_t* jumping_distribution =
                        jumping_distributions[id_right];

                // save old branch lengths
                double d1_old = tree->node_map[id_left ]->d;
                double d2_old = tree->node_map[id_right]->d;
                // generate a proposal
                double eta_old = (d1_old + d2_old)/2.0;
                double  xi_old = (d1_old - d2_old)/2.0;
                double  xi_new = std::min(eta_old, std::max(-eta_old, jumping_distribution->sample(rng, xi_old, 0.01+eta_old)));
                double  d1_new = eta_old + xi_new;
                double  d2_new = eta_old - xi_new;
                print_debug(" d1_old: %f\n",  d1_old);
                print_debug(" d2_old: %f\n",  d2_old);
                print_debug(" d1_new: %f\n",  d1_new);
                print_debug(" d2_new: %f\n",  d2_new);
                print_debug("eta_old: %f\n", eta_old);
                print_debug(" xi_old: %f\n",  xi_old);
                print_debug(" xi_new: %f\n",  xi_new);
                tree->node_map[id_left ]->d = d1_new;
                tree->node_map[id_right]->d = d2_new;

                log_likelihood_new = log_likelihood();

                rho = exp(log_likelihood_new-log_likelihood_ref)
                        *gamma_distribution.pdf(d1_new)/gamma_distribution.pdf(d1_old)
                        *gamma_distribution.pdf(d2_new)/gamma_distribution.pdf(d2_old)
                        *jumping_distribution->p(xi_old, xi_new);
                print_debug("    rho: %f\n",  rho);
                x   = gsl_ran_flat(rng, 0.0, 1.0);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        print_debug("accepted d1 (xi) : %f\n", d1_new);
                        print_debug("accepted d2 (xi) : %f\n", d2_new);
                        update_acceptance(id_right, true);
                        return log_likelihood_new;
                }
                else {
                        // sample rejected
                        tree->node_map[id_left ]->d = d1_old;
                        tree->node_map[id_right]->d = d2_old;
                        print_debug("rejected d1 (xi) : %f\n", d1_new);
                        print_debug("rejected d2 (xi) : %f\n", d2_new);
                        update_acceptance(id_right, false);
                        return log_likelihood_ref;
                }
        }
        double sample(pt_node_t::id_t id, double log_likelihood_ref) {
                log_likelihood_ref = sample_eta(id, log_likelihood_ref);
                log_likelihood_ref = sample_xi (id, log_likelihood_ref);
                print_debug("--------------------------------------\n");
                return log_likelihood_ref;
        }
        void sample(bool print) {
                double log_likelihood_ref = log_likelihood();
                flockfile(stderr);
                std::cerr << "step: "
                          << step
                          << std::endl;
                fflush(stderr);
                funlockfile(stderr);
                // loop over nodes
                for (pt_node_t::id_t id = 0; id < tree->size(); id++) {
                        // if this is not a leaf then sample it
                        if (!tree->node_map[id]->leaf()) {
                                log_likelihood_ref = sample(id, log_likelihood_ref);
                        }
                }
                update_means(print);
                step++;
        }
        void burnin(size_t n, bool print = true) {
                // burn in
                for (size_t i = 0; i < n; i++) {
                        sample(print);
                        update_steps();
                }
                step = 0;
        }
        void sample(size_t n, bool print = true) {
                // sample n times
                for (size_t i = 0; i < n; i++) {
                        sample(print);
                }
        }
        void apply(pt_root_t* tree) {
                // insert means into the tree
                for (pt_root_t::node_map_t::const_iterator is = ++tree->node_map.begin(); is != tree->node_map.end(); is++) {
                        (*is)->d = means[(*is)->id];
                }
        }

        std::vector<double> means;
        std::vector<double> acceptance;
        std::vector<std::vector<double> > history;
private:
        const alignment_t<CODE_TYPE>& alignment;
        exponent_t<CODE_TYPE, ALPHABET_SIZE> alpha;

        gamma_distribution_t gamma_distribution;

        std::vector<jumping_distribution_t*> jumping_distributions;
        double acceptance_rate;

        gsl_rng * rng;
        size_t step;

        pt_root_t* tree;
};

#include <assert.h>
#include <pthread.h>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_pmcmc_hastings_t
{
public:
        pt_pmcmc_hastings_t(size_t n, const pt_metropolis_hastings_t<CODE_TYPE, ALPHABET_SIZE>& mh)
                : means(mh.means) {

                assert(n > 0);

                for (size_t i = 0; i < n; i++) {
                        population.push_back(mh.clone());
                }
        }
        ~pt_pmcmc_hastings_t() {
                for (size_t i = 0; i < population.size(); i++) {
                        delete(population[i]);
                }
        }

        // typedefs
        ////////////////////////////////////////////////////////////////////////////////
        typedef pt_metropolis_hastings_t<CODE_TYPE, ALPHABET_SIZE> pt_sampler_t;

        typedef struct {
                pt_sampler_t* sampler;
                size_t samples;
                size_t burnin;
        } pthread_data_t;

        // sampling methods
        ////////////////////////////////////////////////////////////////////////////////
        static
        void * sample_thread(void* _data) {
                pthread_data_t* data  = (pthread_data_t*)_data;
                pt_sampler_t* sampler = data->sampler;
                const size_t  samples = data->samples;
                const size_t burnin   = data->burnin;

                sampler->burnin(burnin,  false);
                sampler->sample(samples, false);

                return NULL;
        }

        void sample(size_t samples, size_t burnin) {

                pthread_data_t data[population.size()];
                pthread_t threads[population.size()];
                int rc;

                for (size_t i = 0; i < population.size(); i++) {
                        data[i].sampler = population[i];
                        data[i].samples = samples;
                        data[i].burnin  = burnin;
                }

                // sample
                for (size_t i = 0; i < population.size(); i++) {
                        rc = pthread_create(&threads[i], NULL, sample_thread, (void *)&data[i]);
                        if (rc) {
                                std::cerr << "Couldn't create thread." << std::endl;
                                exit(EXIT_FAILURE);
                        }
                }
                // join threads
                for (size_t i = 0; i < population.size(); i++) {
                        rc = pthread_join(threads[i], NULL);
                        if (rc) {
                                std::cerr << "Couldn't join thread." << std::endl;
                                exit(EXIT_FAILURE);
                        }
                }
                update_means();
        }
        void update_means() {
                size_t n = population[0]->history.size();
                size_t t = population[0]->history[0].size();
                for (size_t i = 0; i < t; i++) {
                        for (size_t j = 1; j < n; j++) {
                                double mean = 0.0;
                                for (typename std::vector<pt_sampler_t*>::iterator it = population.begin(); it != population.end(); it++) {
                                        mean += (*it)->history[j][i];
                                }
                                means[j] = mean/(double)population.size();
                        }
                        history.push_back(means);
                }
        }
        void apply(pt_root_t* tree) {
                // insert means into the tree
                for (pt_root_t::node_map_t::const_iterator is = ++tree->node_map.begin(); is != tree->node_map.end(); is++) {
                        (*is)->d = means[(*is)->id];
                }
        }

        std::vector<double> means;
        std::vector<std::vector<double> > history;

private:
        std::vector<pt_sampler_t*> population;
};

#endif /* PHYLOTREE_SAMPLER_HH */
