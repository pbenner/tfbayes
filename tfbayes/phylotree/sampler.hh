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

#ifndef __TFBAYES_PHYLOTREE_SAMPLER_HH__
#define __TFBAYES_PHYLOTREE_SAMPLER_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm> /* std::min */
#include <cmath>
#include <iomanip>
#include <list>
#include <numeric>
#include <vector>

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
#include <tfbayes/phylotree/polynomial.hh>
#include <tfbayes/phylotree/posterior.hh>
#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/polynomial.hh>
#include <tfbayes/utility/progress.hh>
#include <tfbayes/utility/random.hh>
#include <tfbayes/utility/thread-pool.hh>

// misc tools
////////////////////////////////////////////////////////////////////////////////

#define __line_up__  "\x1b[A"
#define __line_del__ "\33[2K\r"

inline
double __rate__(size_t i, size_t n) {
        if (n == 0) return 0.0;
        else        return static_cast<double>(i)/static_cast<double>(n);
}

// proposal distributions
////////////////////////////////////////////////////////////////////////////////
class proposal_distribution_t : public virtual clonable
{
public:
        virtual ~proposal_distribution_t() { }

        proposal_distribution_t* clone() const = 0;

        virtual double p(double d_old, double d_new) const = 0;
        virtual double sample(threaded_rng_t& rng, double d_old) const = 0;
        virtual void increase_jump(double eta) = 0;
        virtual void decrease_jump(double eta) = 0;
};

class normal_proposal_t : public proposal_distribution_t
{
public:
        normal_proposal_t(double sigma_square = 0.2)
                : sigma(std::sqrt(sigma_square)) { }

        normal_proposal_t* clone() const {
                return new normal_proposal_t(*this);
        }

        double p(double d_old, double d_new) const {
                return 1.0;
        }
        double sample(threaded_rng_t& rng, double d_old) const {
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

class gamma_proposal_t : public proposal_distribution_t
{
public:
        gamma_proposal_t(double r, double lambda)
                : gamma_distribution(r, lambda) { }

        gamma_proposal_t* clone() const {
                return new gamma_proposal_t(*this);
        }

        // p(d_new -> d_old)/p(d_old -> d_new) = p(d_old)/p(d_new)
        double p(double d_old, double d_new) const {
                
                return boost::math::pdf(gamma_distribution, d_old)/
                        boost::math::pdf(gamma_distribution, d_new);
        }
        double sample(threaded_rng_t& rng, double d_old) const {
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

// sampling history
////////////////////////////////////////////////////////////////////////////////

struct pt_history_mh_t
{
        typedef std::list<pt_root_t> samples_t;
        typedef std::vector<double> values_t;

        // acceptance rates
        size_t accepted_topologies;
        size_t accepted_lengths;
        size_t steps_topology;
        size_t steps_length;
        // list of tree samples
        samples_t samples;
        // the posterior value for each tree
        values_t values;
};

struct pt_history_mc3_t
{
        pt_history_mc3_t(size_t k)
        : accepted_swaps(k, 0),
          steps         (k, 0)
                { }

        std::vector<size_t> accepted_swaps;
        std::vector<size_t> steps;
};

// the state of a sampler is a phylogenetic tree paired with its
// posterior value
////////////////////////////////////////////////////////////////////////////////
class pt_state_t : std::pair<pt_root_t, double>
{
        typedef std::pair<pt_root_t, double> base_t;
public:
        pt_state_t(const pt_root_t& tree)
                : base_t(tree, 0.0)
                { }
        pt_state_t(const pt_state_t& state)
                : base_t(state)
                { }
        pt_state_t& operator=(const pt_root_t& tree) {
                base_t::first = tree;
                return *this;
        }
        pt_state_t& operator=(double posterior_value) {
                base_t::second = posterior_value;
                return *this;
        }
        const pt_root_t& tree() const {
                return base_t::first;
        }
        pt_root_t& tree() {
                return base_t::first;
        }
        operator pt_root_t&() {
                return base_t::first;
        }
        operator double&() {
                return base_t::second;
        }
};

// phylotree sampler
////////////////////////////////////////////////////////////////////////////////

class pt_sampler_t : public virtual clonable
{
public:
        virtual ~pt_sampler_t() { }

        virtual pt_sampler_t* clone() const = 0;

        // execute sampler
        virtual double operator()(threaded_rng_t& rng, bool verbose = false) = 0;

        // access methods
        virtual const pt_history_mh_t& history() const = 0;
        virtual const pt_state_t& state() const = 0;

        // print status
        virtual void rewind() const = 0;
        virtual void print_progress() const = 0;
};

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class pt_hamiltonian_t : public pt_sampler_t
{
        typedef std::vector<double> vector_t;

        template<typename T>
        struct square
        {
                T operator()(const T& Left, const T& Right) const {   
                        return (Left + Right*Right);
                }
        };
public:
        typedef boost::math::gamma_distribution<> gamma_distribution_t;

        pt_hamiltonian_t(const pt_root_t& tree,
                         const alignment_map_t<AC>& alignment,
                         const std::vector<double>& alpha,
                         const gamma_distribution_t& gamma_distribution,
                         double epsilon,
                         thread_pool_t& thread_pool,
                         double temperature = 1.0)
                : _history_              (), // zero-initialize history
                  _state_                (tree),
                  _thread_pool_          (thread_pool),
                  _alignment_            (alignment),
                  _alpha_                (alpha.begin(), alpha.end()),
                  _gamma_distribution_   (gamma_distribution),
                  _temperature_          (temperature),
                  _epsilon_              (epsilon) {

                // compute the posterior value for the initial tree
                _state_ = log_posterior(_state_);
        }
        pt_hamiltonian_t(const pt_hamiltonian_t& mh)
                : _history_              (mh._history_),
                  _state_                (mh._state_),
                  _thread_pool_          (mh._thread_pool_),
                  _alignment_            (mh._alignment_),
                  _alpha_                (mh._alpha_),
                  _gamma_distribution_   (mh._gamma_distribution_),
                  _temperature_          (mh._temperature_),
                  _epsilon_              (mh._epsilon_) {

        }
        virtual ~pt_hamiltonian_t()
                { }

        virtual pt_hamiltonian_t* clone() const {
                return new pt_hamiltonian_t(*this);
        }

        virtual void rewind() const {
                std::cerr << __line_up__
                          << __line_up__
                          << __line_up__
                          << __line_up__
                          << __line_up__;
        }
        virtual void print_progress() const {
                std::cerr << __line_del__"Sampler with temperature "             << temperature() << std::endl
                          << __line_del__" -> acceptance rate (topology)     : " << __rate__(_history_.accepted_topologies, _history_.steps_topology) << std::endl
                          << __line_del__" -> acceptance rate (branch length): " << __rate__(_history_.accepted_lengths,    _history_.steps_length  ) << std::endl
                          << __line_del__ << std::endl
                          << __line_del__ << std::endl;
        }
        double log_posterior(const pt_root_t& tree) {
                return pt_posterior<AS, AC, PC>(tree, _alignment_, _alpha_, _gamma_distribution_, _thread_pool_);
        }
        pt_marginal_derivative_t log_posterior_derivative(const pt_root_t& tree) {
                return pt_posterior_derivative<AS, AC, PC>(tree, _alignment_, _alpha_, _gamma_distribution_, _thread_pool_);
        }
        void momentum_step(vector_t& p, pt_root_t& q, double epsilon) {
                // current posterior value and derivative
                pt_marginal_derivative_t U = log_posterior_derivative(q);
                for (size_t i = 0; i < p.size(); i++) {
                        p[i] = p[i] - epsilon*U.derivative()[i];
                }
        }
        void position_step(vector_t& p, pt_root_t& q, double epsilon) {
                for (size_t i = 0; i < p.size(); i++) {
                        q[i]->d = std::max(1.0e-20, q[i]->d + epsilon*p[i]);
                }
        }
        void leapfrog(size_t n, vector_t& p, pt_root_t& q) {
                assert(n > 0);
                // half step for momentum
                momentum_step(p, q, _epsilon_/2.0);
                // full steps for momentum and position
                for (size_t i = 0; i < n-1; i++) {
                        position_step(p, q, _epsilon_);
                        momentum_step(p, q, _epsilon_);
                }
                // full step for position
                position_step(p, q, _epsilon_);
                // half step for momentum
                momentum_step(p, q, _epsilon_/2.0);
        }
        void sample_length(threaded_rng_t& rng) {
                size_t n = _state_.tree().n_nodes-1;
                // distributions for drawing random numbers
                boost::random::uniform_01<> uniform;
                // state
                vector_t  p(n);
                pt_root_t q(_state_);
                // current posterior values
                double current_U = _state_;
                double current_K = std::accumulate(p.begin(), p.end(), 0.0, square<double>())/2.0;
                // sample a new momentum
                boost::normal_distribution<> nd(0.0, 1.0);
                for (size_t i = 0; i < n; i++) {
                        p[i] = nd(rng);
                }
                // simulate Hamiltonian dynamics
                leapfrog(1, p, q);
                // Metropolis-Hastings acceptance probability
                double proposed_U = log_posterior(q);
                double proposed_K = std::accumulate(p.begin(), p.end(), 0.0, square<double>())/2.0;
                double rho = std::exp(current_U - proposed_U + current_K - proposed_K);
                if (uniform(rng) < rho) {
                        // accept proposal
                        _state_ = q;
                        _state_ = proposed_U;
                        _history_.accepted_lengths++;
                }
                // update history
                _history_.samples.push_back(_state_);
                _history_.values .push_back(_state_);
                _history_.steps_length++;
        }
        void sample_topology(pt_node_t& node, threaded_rng_t& rng) {
                // can't sample leafs
                if (node.leaf()) return;
                double log_posterior_ref = _state_;
                // distributions for drawing random numbers
                boost::random::bernoulli_distribution<> bernoulli(0.5);
                boost::random::uniform_01<> uniform;
                // this variable indicates the orthant to move to
                size_t which = bernoulli(rng);

                // alter topology
                node.move(which);

                // compute new log likelihood
                const double log_posterior_new = log_posterior(_state_);

                // compute acceptance probability
                const double rho = exp(1.0/_temperature_*(log_posterior_new-log_posterior_ref));
                const double x   = uniform(rng);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        _state_ = log_posterior_new;
                        _history_.accepted_topologies++;
                }
                else {
                        // sample rejected
                        node.move(which);
                }
                _history_.steps_topology++;
        }
        virtual double operator()(threaded_rng_t& rng, bool verbose = false) {
                // sample branch lengths for the full tree
                sample_length(rng);
                // loop over nodes and sample topologies
                for (pt_node_t::nodes_t::const_iterator it = _state_.tree().begin_nodes();
                     it != _state_.tree().end_nodes(); it++) {
                        // skip the root
                        if ((*it)->root()) {
                                continue;
                        }
                        sample_topology(**it, rng);
                }
                // update history
                _history_.samples.push_back(_state_);
                _history_.values .push_back(_state_);
                if (verbose) {
                        print_progress();
                }
                // return posterior value
                return _state_;
        }
        // access methods
        ////////////////////////////////////////////////////////////////////////
        virtual const pt_history_mh_t& history() const {
                return _history_;
        }
        virtual pt_history_mh_t& history() {
                return _history_;
        }
        virtual const pt_state_t& state() const {
                return _state_;
        }
        virtual pt_state_t& state() {
                return _state_;
        }
        const double& temperature() const {
                return _temperature_;
        }
        double& temperature() {
                return _temperature_;
        }
protected:
        // sampler history
        pt_history_mh_t _history_;
        // state of the sampler
        pt_state_t _state_;
        // a thread pool for computing likelihoods
        thread_pool_t& _thread_pool_;
        // alignment data
        const alignment_map_t<AC>& _alignment_;
        // prior distribution and parameters
        exponent_t<AS, PC> _alpha_;
        gamma_distribution_t _gamma_distribution_;
        double _temperature_;
        // integration step size
        double _epsilon_;
};

template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class pt_metropolis_hastings_t : public pt_sampler_t
{
public:
        typedef boost::math::gamma_distribution<> gamma_distribution_t;

        pt_metropolis_hastings_t(const pt_root_t& tree,
                                 const alignment_map_t<AC>& alignment,
                                 const std::vector<double>& alpha,
                                 const gamma_distribution_t& gamma_distribution,
                                 const proposal_distribution_t& proposal_distribution,
                                 thread_pool_t& thread_pool,
                                 double temperature = 1.0)
                : _history_              (), // zero-initialize history
                  _state_                (tree),
                  _thread_pool_          (thread_pool),
                  _alignment_            (alignment),
                  _alpha_                (alpha.begin(), alpha.end()),
                  _gamma_distribution_   (gamma_distribution),
                  _temperature_          (temperature),
                  _proposal_distribution_(proposal_distribution.clone()) {

                // compute the posterior value for the initial tree
                _state_ = log_posterior();
        }
        pt_metropolis_hastings_t(const pt_metropolis_hastings_t& mh)
                : _history_              (mh._history_),
                  _state_                (mh._state_),
                  _thread_pool_          (mh._thread_pool_),
                  _alignment_            (mh._alignment_),
                  _alpha_                (mh._alpha_),
                  _gamma_distribution_   (mh._gamma_distribution_),
                  _temperature_          (mh._temperature_),
                  _proposal_distribution_(mh._proposal_distribution_->clone()) {

        }
        virtual ~pt_metropolis_hastings_t() {
                delete(_proposal_distribution_);
        }

        virtual pt_metropolis_hastings_t* clone() const {
                return new pt_metropolis_hastings_t(*this);
        }

        virtual void rewind() const {
                std::cerr << __line_up__
                          << __line_up__
                          << __line_up__
                          << __line_up__
                          << __line_up__;
        }
        virtual void print_progress() const {
                std::cerr << __line_del__"Sampler with temperature "             << temperature() << std::endl
                          << __line_del__" -> acceptance rate (topology)     : " << __rate__(_history_.accepted_topologies, _history_.steps_topology) << std::endl
                          << __line_del__" -> acceptance rate (branch length): " << __rate__(_history_.accepted_lengths,    _history_.steps_length  ) << std::endl
                          << __line_del__ << std::endl
                          << __line_del__ << std::endl;
        }
        double log_posterior() {
                return pt_posterior<AS, AC, PC>(_state_.tree(), _alignment_, _alpha_, _gamma_distribution_, _thread_pool_);
        }
        void sample_length(pt_node_t& node, threaded_rng_t& rng) {
                double log_posterior_ref = _state_;
                // distributions for drawing random numbers
                boost::random::uniform_01<> uniform;

                // generate a proposal
                double d_old = node.d;
                double d_new = _proposal_distribution_->sample(rng, d_old);
                if (d_new < 0.0) {
                        // the posterior does not have support on
                        // negative branch lengths, so reject immediately
                        return;
                }
                node.d = d_new;

                // compute new log likelihood
                const double log_posterior_new = log_posterior();

                // compute acceptance probability
                const double rho = exp(1.0/_temperature_*(log_posterior_new-log_posterior_ref))
                        *_proposal_distribution_->p(d_old, d_new);
                const double x   = uniform(rng);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        _state_ = log_posterior_new;
                        _history_.accepted_lengths++;
                }
                else {
                        // sample rejected
                        node.d = d_old;
                }
                _history_.steps_length++;
        }
        void sample_topology(pt_node_t& node, threaded_rng_t& rng) {
                // can't sample leafs
                if (node.leaf()) return;
                double log_posterior_ref = _state_;
                // distributions for drawing random numbers
                boost::random::bernoulli_distribution<> bernoulli(0.5);
                boost::random::uniform_01<> uniform;
                // this variable indicates the orthant to move to
                size_t which = bernoulli(rng);

                // alter topology
                node.move(which);

                // compute new log likelihood
                const double log_posterior_new = log_posterior();

                // compute acceptance probability
                const double rho = exp(1.0/_temperature_*(log_posterior_new-log_posterior_ref));
                const double x   = uniform(rng);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        _state_ = log_posterior_new;
                        _history_.accepted_topologies++;
                }
                else {
                        // sample rejected
                        node.move(which);
                }
                _history_.steps_topology++;
        }
        virtual double operator()(threaded_rng_t& rng, bool verbose = false) {
                // loop over nodes
                for (pt_node_t::nodes_t::const_iterator it = _state_.tree().begin_nodes();
                     it != _state_.tree().end_nodes(); it++) {
                        // skip the root
                        if ((*it)->root()) {
                                continue;
                        }
                        sample_length  (**it, rng);
                        sample_topology(**it, rng);
                }
                // update history
                _history_.samples.push_back(_state_);
                _history_.values .push_back(_state_);
                if (verbose) {
                        print_progress();
                }
                // return posterior value
                return _state_;
        }
        // access methods
        ////////////////////////////////////////////////////////////////////////
        virtual const pt_history_mh_t& history() const {
                return _history_;
        }
        virtual pt_history_mh_t& history() {
                return _history_;
        }
        virtual const pt_state_t& state() const {
                return _state_;
        }
        virtual pt_state_t& state() {
                return _state_;
        }
        const double& temperature() const {
                return _temperature_;
        }
        double& temperature() {
                return _temperature_;
        }
protected:
        // sampler history
        pt_history_mh_t _history_;
        // state of the sampler
        pt_state_t _state_;
        // a thread pool for computing likelihoods
        thread_pool_t& _thread_pool_;
        // alignment data
        const alignment_map_t<AC>& _alignment_;
        // prior distribution and parameters
        exponent_t<AS, PC> _alpha_;
        gamma_distribution_t _gamma_distribution_;
        double _temperature_;
        // metropolis-hastings proposal distribution
        proposal_distribution_t* _proposal_distribution_;
};

template <typename T>
class pt_mc3_t : public pt_sampler_t
{
        typedef std::vector<double> vector_t;
public:
        pt_mc3_t(const vector_t& temperatures, const T& mh)
                : _history_     (temperatures.size()),
                  _thread_pool_ (temperatures.size()),
                  uniform_int   (0, temperatures.size()-1) {
                assert(temperatures.size() > 0);
                assert(temperatures[0] == 1.0);
                for (size_t i = 0; i < temperatures.size(); i++) {
                        _population_.push_back(mh.clone());
                        _population_[i].temperature() = temperatures[i];
                }
        }
        pt_mc3_t(const pt_mc3_t& pt_mc3)
                : _history_     (pt_mc3._history_),
                  _thread_pool_ (pt_mc3._thread_pool_),
                  uniform_int   (pt_mc3.uniform_int) {
                for (size_t i = 0; i < pt_mc3._population_.size(); i++) {
                        _population_.push_back(pt_mc3._population_[i].clone());
                }
        }

        virtual pt_mc3_t* clone() const {
                return new pt_mc3_t(*this);
        }

        virtual void rewind() const {
                for (size_t i = 0; i < _population_.size(); i++) {
                        _population_[i].rewind();
                }
                std::cerr << __line_up__
                          << __line_up__;
                for (size_t i = 0; i < _population_.size(); i++) {
                        std::cerr << __line_up__;
                }
        }
        virtual void print_progress() const {
                for (size_t i = 0; i < _population_.size(); i++) {
                        _population_[i].print_progress();
                }
                std::cerr << __line_del__"MC3 swap acceptance rates:" << std::endl;
                for (size_t i = 0; i < _population_.size(); i++) {
                        std::cerr << boost::format(__line_del__" -> sampler %3d: %f\n")
                                % i % __rate__(_history_.accepted_swaps[i], _history_.steps[i]);
                }
                std::cerr << std::endl;
        }
        virtual double operator()(threaded_rng_t& rng, bool verbose = false) {
                using std::swap;
                // future log posterior values
                future_vector_t<double> futures(_population_.size());
                // execute the metropolis algorithm on each chain
                for (size_t i = 0; i < _population_.size(); i++) {
                        boost::function<double ()> f = boost::bind(&pt_sampler_t::operator(),
                                                                   boost::ref(_population_[i]),
                                                                   boost::ref(rng), false);
                        // use local thread pool to execute samplers
                        futures[i] = _thread_pool_.schedule(f);
                }
                // select two chains at random
                const size_t i = uniform_int(rng);
                const size_t j = uniform_int(rng);
                if (i != j) {
                        // update history
                        _history_.steps[i]++;
                        _history_.steps[j]++;
                        // get probabilities
                        const double pi = futures[i].get();
                        const double pj = futures[j].get();
                        const double ti = _population_[i].temperature();
                        const double tj = _population_[j].temperature();
                        // metropolis probability for accepting the swap
                        const double r  = std::min(1.0, std::exp(pi/tj + pj/ti - pi/ti - pj/tj));
                        if (uniform_01(rng) <= r) {
                                // update history
                                _history_.accepted_swaps[i]++;
                                _history_.accepted_swaps[j]++;
                                if (verbose) {
                                        print_progress();
                                }
                                // swap states of the two chains
                                swap(_population_[i].state(),
                                     _population_[j].state());
                        }
                }
                // wait for all processes to finish
                futures.wait();
                // return log posterior
                return futures[0].get();
        }
        // access methods
        ////////////////////////////////////////////////////////////////////////
        virtual const pt_history_mh_t& history() const {
                return _population_[0].history();
        }
        virtual const pt_state_t& state() const {
                return _population_[0].state();
        }
protected:
        // mc3 specific history
        pt_history_mc3_t _history_;
        // a population of samplers
        boost::ptr_vector<T> _population_;
        // a local thread pool
        thread_pool_t _thread_pool_;
        // distribution for randomly selecting two chains
        boost::random::uniform_int_distribution<> uniform_int;
        // distribution for the metropolis update
        boost::random::uniform_01<> uniform_01;
};

class pt_pmcmc_t : public pt_sampler_t
{
public:
        pt_pmcmc_t(size_t n, const pt_sampler_t& mh)
                : _thread_pool_(n) {

                for (size_t i = 0; i < n; i++) {
                        _population_.push_back(mh.clone());
                }
        }
        pt_pmcmc_t(const pt_pmcmc_t& pmcmc)
                : _thread_pool_(pmcmc._thread_pool_) {

                for (size_t i = 0; i < pmcmc._population_.size(); i++) {
                        _population_.push_back(pmcmc._population_[i].clone());
                }
        }

        virtual pt_pmcmc_t* clone() const {
                return new pt_pmcmc_t(*this);
        }

        virtual void rewind() const {
                for (size_t i = 0; i < _population_.size(); i++) {
                        _population_[i].rewind();
                        std::cerr << __line_up__
                                  << __line_up__
                                  << __line_up__;
                }
        }
        virtual void print_progress() const {
                for (size_t i = 0; i < _population_.size(); i++) {
                        std::cerr << std::endl
                                  << "|======================================|"
                                  << std::endl << std::endl;
                        _population_[i].print_progress();
                }
        }
        virtual double operator()(threaded_rng_t& rng, bool verbose = false) {
                // future log posterior values
                future_vector_t<double> futures(_population_.size());
                // loop over population and execute samplers
                for (size_t i = 0; i < _population_.size(); i++) {
                        boost::function<double ()> f = boost::bind(&pt_sampler_t::operator(),
                                                                   boost::ref(_population_[i]),
                                                                   boost::ref(rng), false);
                        // use local thread pool to execute samplers
                        futures[i] = _thread_pool_.schedule(f);
                }
                futures.wait();
                return futures[0].get();
        }
        virtual void operator()(size_t n, threaded_rng_t& rng, bool verbose = false) {
                for (size_t i = 0; i < n; i++) {
                        operator()(rng, verbose);
                        if (verbose && i != 0)
                                rewind();
                        if (verbose) {
                                print_progress();
                                std::cerr << progress_t((i+1.0)/(double)n);
                        }
                }
        }
        // access methods
        ////////////////////////////////////////////////////////////////////////
        virtual const pt_history_mh_t& history() const {
                return _population_[0].history();
        }
        virtual const pt_history_mh_t& history(size_t i) const {
                return _population_[i].history();
        }
        virtual const pt_state_t& state() const {
                return _population_[0].state();
        }
        size_t size() const {
                return _population_.size();
        }
protected:
        // a local thread pool
        thread_pool_t _thread_pool_;
        // a population of samplers
        boost::ptr_vector<pt_sampler_t> _population_;
};

#endif /* __TFBAYES_PHYLOTREE_SAMPLER_HH__ */
