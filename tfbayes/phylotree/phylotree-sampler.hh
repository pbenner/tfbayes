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
#include <limits>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/thread.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/polynomial.hh>
#include <tfbayes/utility/progress.hh>
#include <tfbayes/utility/thread-pool.hh>

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

// data structure optimized for computing likelihoods
template <typename AC = alphabet_code_t>
class alignment_map_t : public std::map<std::vector<AC>, double>
{
        typedef std::map<std::vector<AC>, double> base_t;
public:
        explicit alignment_map_t(const alignment_t<AC>& alignment)
                : base_t() {
                for (typename alignment_t<AC>::const_iterator it = alignment.begin();
                     it != alignment.end(); it++) {
                        base_t::operator[](*it) += 1.0;
                }
        }
};

template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class pt_metropolis_hastings_t : public virtual clonable
{
public:
        typedef std::vector<double> vector_t;
        typedef boost::math::gamma_distribution<> gamma_distribution_t;
        typedef boost::unique_future<double> future_t;
        typedef boost::ptr_vector<future_t> future_vector_t;

        pt_metropolis_hastings_t(const pt_root_t& tree,
                                 const alignment_map_t<AC>& alignment,
                                 const std::vector<double>& alpha,
                                 const gamma_distribution_t& gamma_distribution,
                                 const jumping_distribution_t& jumping_distribution,
                                 thread_pool_t& thread_pool,
                                 double temperature = 1.0)
                : _step_              (0),
                  _tree_              (tree),
                  _thread_pool_       (thread_pool),
                  _alignment_         (alignment),
                  _alpha_             (alpha.begin(), alpha.end()),
                  _gamma_distribution_(gamma_distribution),
                  _temperature_       (temperature) {

                // initialize jumping distributions
                for (pt_node_t::id_t id = 0; id < tree.n_nodes; id++) {
                        jumping_distributions.push_back(jumping_distribution.clone());
                }
        }
        pt_metropolis_hastings_t(const pt_metropolis_hastings_t& mh)
                : _log_posterior_history_(mh._log_posterior_history_),
                  _samples_              (mh._samples_),
                  _step_                 (mh._step_),
                  _tree_                 (mh._tree_),
                  _thread_pool_          (mh._thread_pool_),
                  _alignment_            (mh._alignment_),
                  _alpha_                (mh._alpha_),
                  _gamma_distribution_   (mh._gamma_distribution_),
                  _temperature_          (mh._temperature_),
                  jumping_distributions  (mh.jumping_distributions) {

                // initialize nodes and jumping distributions
                for (pt_node_t::id_t id = 0; id < _tree_.n_nodes; id++) {
                        jumping_distributions[id] = mh.jumping_distributions[id]->clone();
                }
        }
        virtual ~pt_metropolis_hastings_t() {
                // free jumping distributions
                for (pt_node_t::id_t id = 0; id < _tree_.n_nodes; id++) {
                        delete(jumping_distributions[id]);
                }
        }

        virtual pt_metropolis_hastings_t* clone() const {
                return new pt_metropolis_hastings_t(*this);
        }

        double log_likelihood(const std::vector<AC>& column, double n) {
                return n*pt_marginal_likelihood(_tree_, column, _alpha_);
        }
        double log_posterior() {
                double result = 0;
                // results for each column are stored in a vector of
                // futures (with pre allocated capacity)
                future_vector_t futures(_alignment_.size());
                // launch threads to compute the likelihood
                for (typename alignment_map_t<AC>::const_iterator it = _alignment_.begin();
                     it != _alignment_.end(); it++) {
                        boost::function<double ()> f = boost::bind(&pt_metropolis_hastings_t::log_likelihood, this, it->first, it->second);
                        futures.push_back(new future_t());
                        futures.back() = _thread_pool_.schedule(f);
                }
                for (size_t i = 0; i < futures.size(); i++) {
                        result += futures[i].get();
                }
                // prior on branch lengths
                for (pt_node_t::nodes_t::iterator it = _tree_.nodes.begin();
                     it != _tree_.nodes.end(); it++) {
                        // skip the root
                        if ((*it)->root()) {
                                continue;
                        }
                        result += std::log(boost::math::pdf(_gamma_distribution_, (*it)->d));
                }
                return result;
        }
        void update_history(double log_posterior_ref) {
                _samples_.push_back(_tree_);
                _log_posterior_history_.push_back(log_posterior_ref);
                _step_++;
        }
        double sample_branch(pt_node_t& node, double log_posterior_ref,
                             boost::random::mt19937& rng) {
                // distributions for drawing random numbers
                boost::random::bernoulli_distribution<> bernoulli(0.5);
                boost::random::uniform_01<> uniform;
                // in case the topology is sampled, this variable
                // indicates the orthant to move to
                ssize_t which = -1;

                // generate a proposal
                double d_old = node.d;
                double d_new = jumping_distributions[node.id]->sample(rng, d_old);
                if (!node.leaf() && (d_new < 0.0 || bernoulli(rng))) {
                        // propose new topology
                        which = bernoulli(rng);
                        node.move(which);
                }
                d_new  = std::abs(d_new);
                node.d = d_new;

                // compute new log likelihood
                const double log_posterior_new = log_posterior();

                // compute acceptance probability
                const double rho = exp(1.0/_temperature_*(log_posterior_new-log_posterior_ref))
                        *jumping_distributions[node.id]->p(d_old, d_new);
                const double x   = uniform(rng);
                if (x <= std::min(1.0, rho)) {
                        return log_posterior_new;
                }
                else {
                        // sample rejected
                        node.d = d_old;
                        // if topology changed then switch it back
                        node.move(which);

                        return log_posterior_ref;
                }
        }
        void operator()(boost::random::mt19937& rng, bool verbose = false) {
                double log_posterior_ref = log_posterior();
                // loop over nodes
                for (pt_node_t::nodes_t::iterator it = _tree_.nodes.begin();
                     it != _tree_.nodes.end(); it++) {
                        // skip the root
                        if ((*it)->root()) {
                                continue;
                        }
                        // otherwise sample
                        log_posterior_ref = sample_branch(**it, log_posterior_ref, rng);
                }
                update_history(log_posterior_ref);
        }
        void operator()(size_t n, boost::random::mt19937& rng, bool verbose = true) {
                // sample n times
                for (size_t i = 0; i < n; i++) {
                        if (verbose) {
                                std::cerr << progress_t((i+1.0)/(double)n);
                        }
                        operator()(rng);
                }
                if (verbose) std::cerr << std::endl;
        }
        // access methods
        ////////////////////////////////////////////////////////////////////////
        const double& log_posterior() const {
                return _log_posterior_history_.back();
        }
        const vector_t& log_posterior_history() const {
                return _log_posterior_history_;
        }
        const std::list<pt_root_t>& samples() const {
                return _samples_;
        }
        const size_t& step() const {
                return _step_;
        }
        const pt_root_t& tree() const {
                return _tree_;
        }
        pt_root_t& tree() {
                return _tree_;
        }
        const double& temperature() const {
                return _temperature_;
        }
        double& temperature() {
                return _temperature_;
        }
protected:
        // sampler history
        vector_t _log_posterior_history_;
        std::list<pt_root_t> _samples_;
        // state of the sampler
        size_t _step_;
        pt_root_t _tree_;
        // a pool of threads that might be shared among multiple
        // samplers
        thread_pool_t& _thread_pool_;
        // alignment data
        const alignment_map_t<AC>& _alignment_;
        // prior distribution and parameters
        exponent_t<AS, PC> _alpha_;
        gamma_distribution_t _gamma_distribution_;
        double _temperature_;
        // metropolis proposal distribution
        std::vector<jumping_distribution_t*> jumping_distributions;
};

template <typename T>
class pt_mc3_t
{
        typedef std::vector<double> vector_t;
public:
        pt_mc3_t(const vector_t& temperatures, const T& mh)
                : uniform_int(0, temperatures.size()-1) {
                assert(temperatures.size() > 0);
                assert(temperatures[0] == 1.0);
                for (size_t i = 0; i < temperatures.size(); i++) {
                        _population_.push_back(mh.clone());
                        _population_[i]->temperature() = temperatures[i];
                }
        }
        pt_mc3_t(const pt_mc3_t& pt_mc3)
                : uniform_int(pt_mc3.uniform_int) {
                for (size_t i = 0; i < pt_mc3._population_.size(); i++) {
                        _population_.push_back(pt_mc3._population_[i]->clone());
                }
        }
        virtual ~pt_mc3_t() {
                for (size_t i = 0; i < _population_.size(); i++) {
                        delete(_population_[i]);
                }
        }

        virtual pt_mc3_t* clone() const {
                return new pt_mc3_t(*this);
        }

        void operator()(boost::random::mt19937& rng, bool verbose = false) {
                using std::swap;
                // execute the metropolis algorithm on each chain
                for (size_t i = 0; i < _population_.size(); i++) {
                        _population_[i]->operator()(rng);
                }
                // select two chains at random
                const size_t i = uniform_int(rng);
                const size_t j = uniform_int(rng);
                if (i != j) {
                        // get probabilities
                        const double pi = _population_[i]->log_posterior();
                        const double pj = _population_[j]->log_posterior();
                        const double ti = _population_[i]->temperature();
                        const double tj = _population_[j]->temperature();
                        // metropolis probability for accepting the swap
                        const double r  = std::min(1.0, std::exp(pi/tj + pj/ti - pi/ti - pj/tj));
                        if (uniform_01(rng) <= r) {
                                if (verbose) {
                                        std::cerr << boost::format(
                                                "\x1b[A"   // go up one line
                                                "\33[2K\r" // delete line
                                                "MC3 switched states %i and %i.\n")
                                                % i % j;
                                }
                                // swap states of the two chains
                                swap(_population_[i]->tree(),
                                     _population_[j]->tree());
                        }
                }
        }
        void operator()(size_t n, boost::random::mt19937& rng, bool verbose = false) {
                // sample n times
                for (size_t i = 0; i < n; i++) {
                        if (verbose) {
                                std::cerr << progress_t((i+1.0)/(double)n);
                        }
                        operator()(rng);
                }
                if (verbose) std::cerr << std::endl;
        }
        // access methods
        ////////////////////////////////////////////////////////////////////////
        const vector_t& log_posterior_history() const {
                return _population_[0]->log_posterior_history();
        }
        const std::list<pt_root_t>& samples() const {
                return _population_[0]->samples();
        }
protected:
        std::vector<T*> _population_;
        // distribution for randomly selecting two chains
        boost::random::uniform_int_distribution<> uniform_int;
        // distribution for the metropolis update
        boost::random::uniform_01<> uniform_01;
};

template <typename T>
class pt_pmcmc_t
{
public:
        pt_pmcmc_t(size_t n, const T& mh) {

                for (size_t i = 0; i < n; i++) {
                        _population_.push_back(mh.clone());
                }
        }
        pt_pmcmc_t(const pt_pmcmc_t& pmcmc)
                : _samples_              (pmcmc._samples),
                  _log_posterior_history_(pmcmc._log_posterior_history_) {

                for (size_t i = 0; i < pmcmc._population_.size(); i++) {
                        _population_.push_back(pmcmc._population_[i]->clone());
                }
        }
        virtual ~pt_pmcmc_t() {
                for (size_t i = 0; i < _population_.size(); i++) {
                        delete(_population_[i]);
                }
        }

        void update_samples() {
                // fill the list such that the order of samples is preserved
                std::vector<std::list<pt_root_t>::const_iterator> it_vec;
                for (typename std::vector<T*>::const_iterator it = _population_.begin();
                     it != _population_.end(); it++) {
                        it_vec.push_back((*it)->samples().begin());
                }
                while (it_vec[0] != _population_[0]->samples().end())
                {
                        for (size_t i = 0; i < _population_.size(); i++) {
                                _samples_.push_back(*it_vec[i]);
                                // advance iteration for sampler i
                                it_vec[i]++;
                        }
                }
        }
        void update_history() {
                // reset history
                _log_posterior_history_ = std::list<std::vector<double> >();
                // copy history from population
                for (typename std::vector<T*>::const_iterator it = _population_.begin();
                     it != _population_.end(); it++) {
                        _log_posterior_history_.push_back((*it)->log_posterior_history());
                }
        }
        void operator()(boost::random::mt19937& rng, bool verbose = false) {
                for (size_t j = 0; j < _population_.size(); j++) {
                        _population_[j]->operator()(rng, verbose);
                }
        }
        void operator()(size_t n, boost::random::mt19937& rng, bool verbose = false) {
                if (verbose) std::cerr << std::endl << std::endl;
                for (size_t i = 0; i < n; i++) {
                        operator()(rng, verbose);
                        if (verbose) {
                                std::cerr << progress_t((i+1.0)/(double)n);
                        }
                }
                if (verbose) std::cerr << std::endl;
                update_samples();
                update_history();
        }
        // access methods
        ////////////////////////////////////////////////////////////////////////
        const std::list<std::vector<double> >& log_posterior_history() const {
                return _log_posterior_history_;
        }
        const std::list<pt_root_t>& samples() const {
                return _samples_;
        }
protected:
        std::list<pt_root_t> _samples_;
        std::list<std::vector<double> > _log_posterior_history_;
        std::vector<T*> _population_;
};

#endif /* __TFBAYES_PHYLOTREE_PHYLOTREE_SAMPLER_HH__ */
