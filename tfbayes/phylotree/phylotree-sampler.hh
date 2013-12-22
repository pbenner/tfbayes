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

#ifndef __TFBAYES_PHYLOTREE_PHYLOTREE_SAMPLER_HH__
#define __TFBAYES_PHYLOTREE_PHYLOTREE_SAMPLER_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <vector>
#include <cmath>
#include <algorithm> /* std::min */
#include <iomanip>
#include <sys/time.h>
#include <limits>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/polynomial.hh>
#include <tfbayes/utility/progress.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

class jumping_distribution_t : public virtual clonable
{
public:
        virtual ~jumping_distribution_t() { }

        jumping_distribution_t* clone() const = 0;

        virtual double p(double d_old, double d_new) const = 0;
        virtual double sample(boost::random::mt19937& rng, double d_old) const = 0;
        virtual void increase_jump(double eta) = 0;
        virtual void decrease_jump(double eta) = 0;
};

class normal_jump_t : public jumping_distribution_t
{
public:
        normal_jump_t(double sigma_square = 0.2)
                : sigma(std::sqrt(sigma_square)) { }

        normal_jump_t* clone() const {
                return new normal_jump_t(*this);
        }

        double p(double d_old, double d_new) const {
                return 1.0;
        }
        double sample(boost::random::mt19937& rng, double d_old) const {
                boost::random::normal_distribution<> dist(0.0, sigma);
                return d_old+dist(rng);
        }
        void increase_jump(double eta) {
                sigma = sigma+eta;
        }
        void decrease_jump(double eta) {
                sigma = std::max(sigma-eta, 0.0);
        }
private:
        double sigma;
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
                
                return boost::math::pdf(gamma_distribution, d_old)/
                        boost::math::pdf(gamma_distribution, d_new);
        }
        double sample(boost::random::mt19937& rng, double d_old) const {
                boost::random::gamma_distribution<> dist(
                        gamma_distribution.shape(),
                        gamma_distribution.scale());
                return dist(rng);
        }
        void increase_jump(double eta) {
        }
        void decrease_jump(double eta) {
        }
private:
        boost::math::gamma_distribution<> gamma_distribution;
};

template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class pt_metropolis_hastings_t : public virtual clonable
{
public:
        typedef std::map<std::vector<AC>, double> alignment_map_t;
        typedef boost::math::gamma_distribution<> gamma_distribution_t;

        pt_metropolis_hastings_t(const pt_root_t& tree,
                                 const alignment_t<AC>& alignment,
                                 const exponent_t<AS, PC>& alpha,
                                 const gamma_distribution_t& gamma_distribution,
                                 const jumping_distribution_t& jumping_distribution)
                : samples(),
                  alpha(alpha),
                  gamma_distribution(gamma_distribution),
                  step(0),
                  tree(tree) {

                // sleep for a millisecond to make sure that we get
                // a unique seed
                boost::this_thread::sleep(boost::posix_time::milliseconds(1));
                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                rng.seed(tv.tv_sec*tv.tv_usec);

                // initialize jumping distributions
                for (pt_node_t::id_t id = 0; id < tree.n_nodes; id++) {
                        jumping_distributions.push_back(jumping_distribution.clone());
                }
                // fill alignment map
                for (typename alignment_t<AC>::const_iterator it = alignment.begin();
                     it != alignment.end(); it++) {
                        this->alignment[*it] += 1.0;
                }
        }
        pt_metropolis_hastings_t(const pt_metropolis_hastings_t& mh)
                : samples(mh.samples),
                  alignment(mh.alignment),
                  alpha(mh.alpha),
                  gamma_distribution(mh.gamma_distribution),
                  jumping_distributions(mh.jumping_distributions),
                  step(mh.step),
                  tree(mh.tree)
                {

                // sleep for a millisecond to make sure that we get
                // a unique seed
                boost::this_thread::sleep(boost::posix_time::milliseconds(1));
                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                rng.seed(tv.tv_sec*tv.tv_usec);

                // initialize nodes and jumping distributions
                for (pt_node_t::id_t id = 0; id < tree.n_nodes; id++) {
                        jumping_distributions[id] = mh.jumping_distributions[id]->clone();
                }
        }
        virtual ~pt_metropolis_hastings_t() {
                // free jumping distributions
                for (pt_node_t::id_t id = 0; id < tree.n_nodes; id++) {
                        delete(jumping_distributions[id]);
                }
        }

        virtual pt_metropolis_hastings_t* clone() const {
                return new pt_metropolis_hastings_t(*this);
        }

        double log_posterior() {
                double result = 0;
                // likelihood
                for (typename alignment_map_t::const_iterator it = alignment.begin(); it != alignment.end(); it++) {
                        result += it->second*pt_marginal_likelihood(tree, it->first, alpha);
                }
                // prior on branch lengths
                for (pt_node_t::nodes_t::iterator it = tree.nodes.begin();
                     it != tree.nodes.end(); it++) {
                        // skip the root
                        if ((*it)->root()) {
                                continue;
                        }
                        result += std::log(boost::math::pdf(gamma_distribution, (*it)->d));
                }
                return result;
        }
        void update_samples() {
                samples.push_back(tree);
        }
        double sample_branch(pt_node_t& node, double log_posterior_ref) {
                // random number generators
                boost::random::bernoulli_distribution<> bernoulli(0.5);
                boost::random::uniform_01<> uniform;
                double rho;
                double x;
                double log_posterior_new;
                int which = -1;

                // generate a proposal
                double d_old = node.d;
                double d_new = jumping_distributions[node.id]->sample(rng, d_old);
                if (!node.leaf() && (d_new < 0.0 || bernoulli(rng))) {
                        // propose new topology
                        which = bernoulli(rng);
                        switch (which) {
                        case 0: node.move_a(); break;
                        case 1: node.move_b(); break;
                        default: break;
                        }
                }
                d_new  = std::abs(d_new);
                node.d = d_new;

                // compute new log likelihood
                log_posterior_new = log_posterior();

                // compute acceptance probability
                rho = exp(log_posterior_new-log_posterior_ref)
                        *jumping_distributions[node.id]->p(d_old, d_new);
                x   = uniform(rng);
                if (x <= std::min(1.0, rho)) {
                        return log_posterior_new;
                }
                else {
                        // sample rejected
                        node.d = d_old;
                        // if topology changed then switch it back
                        switch (which) {
                        case 0: node.move_a(); break;
                        case 1: node.move_b(); break;
                        default: break;
                        }
                        return log_posterior_ref;
                }
        }
        void print_progress(size_t i, size_t n) {
                flockfile(stderr);
                fflush(stderr);
                std::cerr << progress_t((i+1.0)/(double)n);
                funlockfile(stderr);
        }
        virtual void operator()() {
                double log_posterior_ref = log_posterior();
                // loop over nodes
                for (pt_node_t::nodes_t::iterator it = tree.nodes.begin();
                     it != tree.nodes.end(); it++) {
                        // skip the root
                        if ((*it)->root()) {
                                continue;
                        }
                        // otherwise sample
                        log_posterior_ref = sample_branch(**it, log_posterior_ref);
                }
                update_samples();
                log_posterior_history.push_back(log_posterior_ref);
                step++;
        }
        virtual void operator()(size_t n, bool progress = true) {
                // sample n times
                for (size_t i = 0; i < n; i++) {
                        if (progress) print_progress(i, n);
                        operator()();
                }
                if (progress) std::cerr << std::endl;
        }

        std::vector<double> log_posterior_history;
        std::list<pt_root_t> samples;
protected:
        alignment_map_t alignment;
        exponent_t<AS, PC> alpha;

        gamma_distribution_t gamma_distribution;

        std::vector<jumping_distribution_t*> jumping_distributions;

        boost::random::mt19937 rng;
        size_t step;

        pt_root_t tree;
};

template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class pt_pmcmc_hastings_t
{
public:
        pt_pmcmc_hastings_t(size_t n, const pt_metropolis_hastings_t<AS, AC, PC>& mh) {

                for (size_t i = 0; i < n; i++) {
                        population.push_back(mh.clone());
                }
        }
        pt_pmcmc_hastings_t(const pt_pmcmc_hastings_t& pmcmc) {
                for (size_t i = 0; i < pmcmc.population.size(); i++) {
                        population.push_back(pmcmc.population[i]->clone());
                }
                // free tree samples
                for (std::list<pt_root_t>::const_iterator it = pmcmc.samples.begin();
                     it != pmcmc.samples.end(); it++) {
                        samples.push_back(*it);
                }
        }
        ~pt_pmcmc_hastings_t() {
                for (size_t i = 0; i < population.size(); i++) {
                        delete(population[i]);
                }
        }

        // typedefs
        ////////////////////////////////////////////////////////////////////////////////
        typedef pt_metropolis_hastings_t<AS, AC, PC> pt_sampler_t;

        // sampling methods
        ////////////////////////////////////////////////////////////////////////////////
        void operator()(size_t samples, bool verbose = true) {

                std::vector<boost::thread*> threads(population.size());

                // sample
                for (size_t i = 0; i < population.size(); i++) {
                        threads[i] = new boost::thread(boost::ref(*population[i]), samples, i==0 && verbose);
                }
                // join threads
                for (size_t i = 0; i < population.size(); i++) {
                        threads[i]->join();
                }
                for (size_t i = 0; i < population.size(); i++) {
                        delete(threads[i]);
                }
                update_samples();
                update_history();
        }
        void update_samples() {
                // fill the list such that the order of samples is preserved
                std::vector<std::list<pt_root_t>::const_iterator> it_vec;
                for (typename std::vector<pt_sampler_t*>::const_iterator it = population.begin();
                     it != population.end(); it++) {
                        it_vec.push_back((*it)->samples.begin());
                }
                while (it_vec[0] != population[0]->samples.end())
                {
                        for (size_t i = 0; i < population.size(); i++) {
                                samples.push_back(*it_vec[i]);
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

#endif /* __TFBAYES_PHYLOTREE_PHYLOTREE_SAMPLER_HH__ */
