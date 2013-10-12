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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <vector>
#include <cmath>
#include <algorithm> /* std::min */
#include <iomanip>
#include <sys/time.h>

#include <gsl/gsl_randist.h>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/exception/exception.h>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/utility/distribution.hh>
#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/polynomial.hh>
#include <tfbayes/utility/progress.hh>

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
        pt_metropolis_hastings_t(const pt_root_t& tree,
                                 const alignment_t<CODE_TYPE>& alignment,
                                 const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
                                 double r, double lambda,
                                 const jumping_distribution_t& jumping_distribution,
                                 double acceptance_rate = 0.7)
                : acceptance(tree->n_nodes, 0.0),
                  samples(),
                  alignment(alignment),
                  alpha(alpha),
                  gamma_distribution(r, lambda),
                  acceptance_rate(acceptance_rate),
                  step(0),
                  tree(tree) {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize jumping distributions
                for (pt_node_t::id_t id = 0; id < tree->n_nodes; id++) {
                        jumping_distributions.push_back(jumping_distribution.clone());
                }
        }
        pt_metropolis_hastings_t(const pt_metropolis_hastings_t& mh)
                : acceptance(mh.acceptance),
                  samples(),
                  alignment(mh.alignment),
                  alpha(mh.alpha),
                  gamma_distribution(mh.gamma_distribution),
                  jumping_distributions(mh.jumping_distributions),
                  acceptance_rate(mh.acceptance_rate),
                  step(mh.step),
                  tree(mh.tree)
                {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize nodes and jumping distributions
                for (pt_node_t::id_t id = 0; id < tree->n_nodes; id++) {
                        jumping_distributions[id] = mh.jumping_distributions[id]->clone();
                }
                // clone tree samples
                for (std::list<pt_root_t*>::const_iterator it = mh.samples.begin();
                     it != mh.samples.end(); it++) {
                        samples.push_back((*it)->clone());
                }
        }
        virtual ~pt_metropolis_hastings_t() {
                // free random generator
                gsl_rng_free(rng);
                // free jumping distributions
                for (pt_node_t::id_t id = 0; id < tree->n_nodes; id++) {
                        delete(jumping_distributions[id]);
                }
        }

        virtual pt_metropolis_hastings_t* clone() const {
                return new pt_metropolis_hastings_t(*this);
        }

        double log_posterior(const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& polynomial) {
                // loop over monomials
                double result = -HUGE_VAL;
                for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator ut = polynomial.begin();
                     ut != polynomial.end(); ut++) {
                        result = logadd(result, log(ut->coefficient()) + mbeta_log(ut->exponent(), alpha) - mbeta_log(alpha));
                }
                return result;
        }
        double log_posterior() {
                double result = 0;
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        const polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial = pt_polynomial<CODE_TYPE, ALPHABET_SIZE>(tree, *it);
                        result += log_posterior(polynomial);
                }
                return result;
        }
        void update_samples() {
                samples.push_back(tree->clone());
        }
        void update_steps() {
                for (pt_node_t::id_t id = 1; id < tree->n_nodes; id++) {
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
        double sample_branch(pt_node_t* node, double log_posterior_ref) {
                double rho;
                double x;
                double log_posterior_new;
                int which = -1;

                // generate a proposal
                double d_old = node->d;
                double d_new = jumping_distributions[node->id]->sample(rng, d_old);
                if (!node->leaf() && (d_new < 0.0 ||  gsl_ran_bernoulli(rng, 0.5))) {
                        print_debug("proposing new topology\n");
                        // propose new topology
                        which = gsl_ran_bernoulli(rng, 0.5);
                        switch (which) {
                        case 0: node.move_a(); break;
                        case 1: node.move_b(); break;
                        default: break;
                        }
                }
                if (d_new < 0.0) {
                        d_new = -d_new;
                }
                node->d = d_new;

                // compute new log likelihood
                log_posterior_new = log_posterior();

                print_debug("proposal for node %d: %f\n", (int)node->id, d_new);
                print_debug("likelihood reference: %f\n", log_posterior_ref);
                print_debug("likelihood proposal : %f\n", log_posterior_new);

                // compute acceptance probability
                rho = exp(log_posterior_new-log_posterior_ref)
                        *gamma_distribution.pdf(d_new)/gamma_distribution.pdf(d_old)
                        *jumping_distributions[node->id]->p(d_old, d_new);
                x   = gsl_ran_flat(rng, 0.0, 1.0);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        print_debug("accepted: %f\n", d_new);
                        update_acceptance(node->id, true);
                        return log_posterior_new;
                }
                else {
                        // sample rejected
                        node->d = d_old;
                        // if topology changed then switch it back
                        switch (which) {
                        case 0: node.move_a(); break;
                        case 1: node.move_b(); break;
                        default: break;
                        }
                        print_debug("rejected: %f\n", d_new);
                        update_acceptance(node.id, false);
                        return log_posterior_ref;
                }
        }
        virtual void generate_sample() {
                double log_posterior_ref = log_posterior();
                // loop over nodes
                for (pt_node_t::nodes_t::iterator it = tree.nodes.begin();
                     it != tree->nodes.end(); it++) {
                        // skip the root
                        if ((*it)->root()) {
                                continue;
                        }
                        // otherwise sample
                        log_posterior_ref = sample_branch(*it, log_posterior_ref);
                }
                update_samples();
                log_posterior_history.push_back(log_posterior_ref);
                step++;
        }
        void print_progress(size_t i, size_t n) {
                flockfile(stderr);
                fflush(stderr);
                std::cerr << progress_t((i+1.0)/(double)n);
                funlockfile(stderr);
        }
        void burnin(size_t n, bool progress = true) {
                // burn in
                for (size_t i = 0; i < n; i++) {
                        if (progress) print_progress(i, n);
                        generate_sample();
                        update_steps();
                }
                step = 0;
        }
        void sample(size_t n, bool progress = true) {
                // sample n times
                for (size_t i = 0; i < n; i++) {
                        if (progress) print_progress(i, n);
                        generate_sample();
                }
                if (progress) std::cerr << std::endl;
        }

        std::vector<double> acceptance;
        std::vector<double> log_posterior_history;
        std::list<pt_root_t> samples;
protected:
        const alignment_t<CODE_TYPE>& alignment;
        exponent_t<CODE_TYPE, ALPHABET_SIZE> alpha;

        gamma_distribution_t gamma_distribution;

        std::vector<jumping_distribution_t*> jumping_distributions;
        double acceptance_rate;

        gsl_rng * rng;
        size_t step;

        pt_root_t tree;
};

#include <assert.h>
#include <pthread.h>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_pmcmc_hastings_t
{
public:
        pt_pmcmc_hastings_t(size_t n, const pt_metropolis_hastings_t<CODE_TYPE, ALPHABET_SIZE>& mh) {

                assert(n > 0);

                for (size_t i = 0; i < n; i++) {
                        population.push_back(mh.clone());
                }
        }
        pt_pmcmc_hastings_t(const pt_pmcmc_hastings_t& pmcmc) {
                for (size_t i = 0; i < pmcmc.population.size(); i++) {
                        population.push_back(pmcmc.population[i]->clone());
                }
                // free tree samples
                for (std::list<pt_root_t*>::const_iterator it = pmcmc.samples.begin();
                     it != pmcmc.samples.end(); it++) {
                        samples.push_back((*it)->clone());
                }
        }
        ~pt_pmcmc_hastings_t() {
                for (size_t i = 0; i < population.size(); i++) {
                        delete(population[i]);
                }
                // free tree samples
                for (std::list<pt_root_t*>::const_iterator it = samples.begin();
                     it != samples.end(); it++) {
                        (*it)->destroy();
                }
        }

        // typedefs
        ////////////////////////////////////////////////////////////////////////////////
        typedef pt_metropolis_hastings_t<CODE_TYPE, ALPHABET_SIZE> pt_sampler_t;

        typedef struct {
                pt_sampler_t* sampler;
                size_t samples;
                size_t burnin;
                bool print_status;
        } pthread_data_t;

        // sampling methods
        ////////////////////////////////////////////////////////////////////////////////
        static
        void * sample_thread(void* _data) {
                pthread_data_t* data  = (pthread_data_t*)_data;
                pt_sampler_t* sampler = data->sampler;
                const size_t  samples = data->samples;
                const size_t burnin   = data->burnin;

                sampler->burnin(burnin,  data->print_status);
                sampler->sample(samples, data->print_status);

                return NULL;
        }
        void sample(size_t samples, size_t burnin) {

                pthread_data_t* data = new pthread_data_t[population.size()];
                pthread_t threads[population.size()];
                int rc;

                for (size_t i = 0; i < population.size(); i++) {
                        data[i].sampler      = population[i];
                        data[i].samples      = samples;
                        data[i].burnin       = burnin;
                        data[i].print_status = (i == 0);
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
                update_samples();
                update_history();
        }
        void delete_samples() {
                for (typename std::list<pt_root_t>::const_iterator is = samples.begin();
                     is != samples.end(); is++) {
                        (*is)->destroy();
                }
                samples = std::list<pt_root_t*>();
        }
        void update_samples() {
                // first remove all samples from previous runs
                delete_samples();
                // fill the list such that the order of samples is preserved
                std::vector<std::list<pt_root_t*>::const_iterator> it_vec;
                for (typename std::vector<pt_sampler_t*>::const_iterator it = population.begin();
                     it != population.end(); it++) {
                        it_vec.push_back((*it)->samples.begin());
                }
                while (it_vec[0] != population[0]->samples.end())
                {
                        for (size_t i = 0; i < population.size(); i++) {
                                samples.push_back((*it_vec[i])->clone());
                                // advance iteration for sampler i
                                it_vec[i]++;
                        }
                }
        }
        void delete_history() {
                log_posterior_history = std::list<std::vector<double> >();
        }
        void update_history() {
                delete_history();
                for (typename std::vector<pt_sampler_t*>::const_iterator it = population.begin();
                     it != population.end(); it++) {
                        log_posterior_history.push_back((*it)->log_posterior_history);
                }
        }

        std::list<pt_root_t> samples;
        std::list<std::vector<double> > log_posterior_history;
private:
        std::vector<pt_sampler_t*> population;
};

#endif /* PHYLOTREE_SAMPLER_HH */
