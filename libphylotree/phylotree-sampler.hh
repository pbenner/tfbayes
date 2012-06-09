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
//                return std::max(d_old+gsl_ran_gaussian(rng, sigma_square), 0.0);
                return d_old+gsl_ran_gaussian(rng, sigma_square);
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
                  nodes(tree->size()) {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize jumping distributions
                for (pt_node_t::id_t id = 0; id < nodes; id++) {
                        jumping_distributions.push_back(jumping_distribution.clone());
                }

                // initialize trees
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        // add tree
                        pt_root_t* _tree = tree->clone();
                        it.apply(_tree);
                        trees.push_back(_tree);
                        // add polynomial
                        pt_cached_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(_tree);
                        polynomials    .push_back(polynomial);
                        polynomials_tmp.push_back(polynomial);
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
                  nodes(mh.nodes)
                {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize nodes and jumping distributions
                for (pt_node_t::id_t id = 0; id < nodes; id++) {
                        jumping_distributions[id] = mh.jumping_distributions[id]->clone();
                }
                // initialize trees
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        // add tree
                        pt_root_t* _tree = mh.trees[0]->clone();
                        it.apply(_tree);
                        trees.push_back(_tree);
                        // add polynomial
                        pt_cached_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(_tree);
                        polynomials    .push_back(polynomial);
                        polynomials_tmp.push_back(polynomial);
                }
        }
        ~pt_metropolis_hastings_t() {
                // free random generator
                gsl_rng_free(rng);
                // free jumping distributions
                for (pt_node_t::id_t id = 0; id < nodes; id++) {
                        delete(jumping_distributions[id]);
                }
                // free trees
                for (trees_t::iterator it = trees.begin(); it != trees.end(); it++) {
                        (*it)->destroy();
                }
        }

        pt_metropolis_hastings_t* clone() const {
                return new pt_metropolis_hastings_t(*this);
        }

        double log_likelihood() {
                double result = 0;
                for (typename polynomials_t::iterator it = polynomials.begin(); it != polynomials.end(); it++) {
                        // loop over monomials
                        double tmp = -HUGE_VAL;
                        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator ut = it->begin();
                             ut != it->end(); ut++) {
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
                for (pt_node_t::id_t id = 1; id < nodes; id++) {
                        pt_node_t* node = trees[0]->node_map[id];
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
                for (pt_node_t::id_t id = 1; id < nodes; id++) {
                        if (acceptance[id] > acceptance_rate) {
                                jumping_distributions[id]->increase_jump(fabs(acceptance[id]-acceptance_rate));
                        }
                        else {
                                jumping_distributions[id]->decrease_jump(fabs(acceptance[id]-acceptance_rate));
                        }
                }
        }
        void propose(pt_node_t::id_t id, double d_old, double d_new) {

                print_debug("old sample: %f\n", d_old);

                // compute new log likelihood
                for (trees_t::iterator it = trees.begin(); it != trees.end(); it++) {
                        (*it)->node_map[id]->d = d_new;
                }
                for (typename polynomials_t::iterator it = polynomials.begin(); it != polynomials.end(); it++) {
                        it->update_length(id);
                }
        }
        void accept() {
                for (size_t i = 0; i < polynomials.size(); i++) {
                        polynomials_tmp[i] = polynomials[i];
                }
        }
        void reject(pt_node_t::id_t id, double d_old) {
                for (size_t i = 0; i < polynomials.size(); i++) {
                        polynomials[i] = polynomials_tmp[i];
                        trees[i]->node_map[id]->d = d_old;
                }
        }
        std::pair<double,bool> sample(pt_node_t::id_t id, double log_likelihood_ref) {
                double rho;
                double x;
                double log_likelihood_new;

                // generate a proposal
                double d_old = trees[0]->node_map[id]->d;
                double d_new = jumping_distributions[id]->sample(rng, d_old);
                propose(id, d_old, d_new);

                log_likelihood_new = log_likelihood();

                // compute acceptance probability
                rho = exp(log_likelihood_new-log_likelihood_ref)
                        *gamma_distribution.pdf(d_new)/gamma_distribution.pdf(d_old)
                        *jumping_distributions[id]->p(d_old, d_new);
                x   = gsl_ran_flat(rng, 0.0, 1.0);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        accept();
                        print_debug("accepted: %f\n", d_new);
                        return std::pair<double, bool>(log_likelihood_new, true);
                }
                else {
                        // sample rejected
                        reject(id, d_old);
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
                for (pt_node_t::id_t id = 1; id < nodes; id++) {
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

        // for each position of the alignment keep a tree that
        // represents this position and two vectors of polynomials for
        // the likelihood
        typedef std::vector<pt_root_t*> trees_t;
        typedef std::vector<pt_cached_polynomial_t<CODE_TYPE, ALPHABET_SIZE> > polynomials_t;

        pt_node_t::id_t nodes;
        trees_t trees;
        polynomials_t polynomials;
        polynomials_t polynomials_tmp;
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
                  nodes(tree->size()) {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize jumping distributions
                for (pt_node_t::id_t id = 0; id < nodes; id++) {
                        jumping_distributions.push_back(jumping_distribution.clone());
                }

                // initialize trees
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        // add tree
                        pt_root_t* _tree = tree->clone();
                        it.apply(_tree);
                        trees.push_back(_tree);
                        // add polynomial
                        pt_cached_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(_tree);
                        polynomials    .push_back(polynomial);
                        polynomials_tmp.push_back(polynomial);
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
                  nodes(mh.nodes)
                {

                // initialize random generator
                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                rng = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(rng, seed);

                // initialize nodes and jumping distributions
                for (pt_node_t::id_t id = 0; id < nodes; id++) {
                        jumping_distributions[id] = mh.jumping_distributions[id]->clone();
                }
                // initialize trees
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        // add tree
                        pt_root_t* _tree = mh.trees[0]->clone();
                        it.apply(_tree);
                        trees.push_back(_tree);
                        // add polynomial
                        pt_cached_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(_tree);
                        polynomials    .push_back(polynomial);
                        polynomials_tmp.push_back(polynomial);
                }
        }
        ~pt_geometric_hastings_t() {
                // free random generator
                gsl_rng_free(rng);
                // free jumping distributions
                for (pt_node_t::id_t id = 0; id < nodes; id++) {
                        delete(jumping_distributions[id]);
                }
                // free trees
                for (trees_t::iterator it = trees.begin(); it != trees.end(); it++) {
                        (*it)->destroy();
                }
        }

        pt_geometric_hastings_t* clone() const {
                return new pt_geometric_hastings_t(*this);
        }

        double log_likelihood() {
                double result = 0;
                for (typename polynomials_t::iterator it = polynomials.begin(); it != polynomials.end(); it++) {
                        // loop over monomials
                        double tmp = -HUGE_VAL;
                        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator ut = it->begin();
                             ut != it->end(); ut++) {
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
                for (pt_node_t::id_t id = 1; id < nodes; id++) {
                        pt_node_t* node = trees[0]->node_map[id];
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
                for (pt_node_t::id_t id = 1; id < nodes; id++) {
                        if (acceptance[id] > acceptance_rate) {
                                jumping_distributions[id]->increase_jump(fabs(acceptance[id]-acceptance_rate));
                        }
                        else {
                                jumping_distributions[id]->decrease_jump(fabs(acceptance[id]-acceptance_rate));
                        }
                }
        }
        void propose(pt_node_t::id_t id1, pt_node_t::id_t id2, double d1_new, double d2_new) {
                // compute new log likelihood
                for (trees_t::iterator it = trees.begin(); it != trees.end(); it++) {
                        (*it)->node_map[id1]->d = d1_new;
                        (*it)->node_map[id2]->d = d2_new;
                }
                for (typename polynomials_t::iterator it = polynomials.begin(); it != polynomials.end(); it++) {
                        it->update_length(id1, id2);
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
        void accept() {
                for (size_t i = 0; i < polynomials.size(); i++) {
                        polynomials_tmp[i] = polynomials[i];
                }
        }
        void reject(pt_node_t::id_t id1, pt_node_t::id_t id2, double d1_old, double d2_old) {
                for (size_t i = 0; i < polynomials.size(); i++) {
                        polynomials[i] = polynomials_tmp[i];
                        trees[i]->node_map[id1]->d = d1_old;
                        trees[i]->node_map[id2]->d = d2_old;
                }
        }
        double sample_eta(pt_node_t::id_t id, double log_likelihood_ref) {
                double rho;
                double x;
                double log_likelihood_new;

                pt_node_t::id_t id_left  = trees[0]->node_map[id]->left ->id;
                pt_node_t::id_t id_right = trees[0]->node_map[id]->right->id;
                jumping_distribution_t* jumping_distribution =
                        jumping_distributions[id_left];

                // save old branch lengths
                double d1_old = trees[0]->node_map[id]->left ->d;
                double d2_old = trees[0]->node_map[id]->right->d;
                // generate a proposal
                double eta_old = (d1_old + d2_old)/2.0;
                double  xi_old = (d1_old - d2_old)/2.0;
                double eta_new = jumping_distribution->sample(rng, eta_old);
                double  d1_new = eta_new + xi_old;
                double  d2_new = eta_new - xi_old;
                print_debug(" d1_old: %f\n",  d1_old);
                print_debug(" d2_old: %f\n",  d2_old);
                print_debug(" d1_new: %f\n",  d1_new);
                print_debug(" d2_new: %f\n",  d2_new);
                print_debug("eta_old: %f\n", eta_old);
                print_debug(" xi_old: %f\n",  xi_old);
                print_debug("eta_new: %f\n", eta_new);
                propose(id_left, id_right, d1_new, d2_new);

                log_likelihood_new = log_likelihood();

                // compute acceptance probability
                if (d1_new < 0.0 || d2_new < 0.0) {
                        rho = 0.0;
                }
                else {
                        rho = exp(log_likelihood_new-log_likelihood_ref)
                                *gamma_distribution.pdf(d1_new)/gamma_distribution.pdf(d1_old)
                                *gamma_distribution.pdf(d2_new)/gamma_distribution.pdf(d2_old)
                                *jumping_distribution->p(eta_old, eta_new);
                }
                x = gsl_ran_flat(rng, 0.0, 1.0);
                print_debug("rho: %f\n", rho);
                print_debug("rho: %f\n", std::min(1.0, rho));
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        accept();
                        print_debug("accepted d1 (eta): %f\n", d1_new);
                        print_debug("accepted d2 (eta): %f\n", d2_new);
                        update_acceptance(id_left, true);
                        return log_likelihood_new;
                }
                else {
                        // sample rejected
                        reject(id_left, id_right, d1_old, d2_old);
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

                pt_node_t::id_t id_left  = trees[0]->node_map[id]->left ->id;
                pt_node_t::id_t id_right = trees[0]->node_map[id]->right->id;
                jumping_distribution_t* jumping_distribution =
                        jumping_distributions[id_right];

                // save old branch lengths
                double d1_old = trees[0]->node_map[id]->left ->d;
                double d2_old = trees[0]->node_map[id]->right->d;
                // generate a proposal
                double eta_old = (d1_old + d2_old)/2.0;
                double  xi_old = (d1_old - d2_old)/2.0;
                double  xi_new = jumping_distribution->sample(rng, xi_old);
                double  d1_new = std::max(0.0, eta_old + xi_new);
                double  d2_new = std::max(0.0, eta_old - xi_new);
                print_debug(" d1_old: %f\n",  d1_old);
                print_debug(" d2_old: %f\n",  d2_old);
                print_debug(" d1_new: %f\n",  d1_new);
                print_debug(" d2_new: %f\n",  d2_new);
                print_debug("eta_old: %f\n", eta_old);
                print_debug(" xi_old: %f\n",  xi_old);
                print_debug(" xi_new: %f\n",  xi_new);
                propose(id_left, id_right, d1_new, d2_new);

                log_likelihood_new = log_likelihood();

                // compute acceptance probability
                if (d1_new < 0.0 || d2_new < 0.0) {
                        rho = 0.0;
                }
                else {
                        rho = exp(log_likelihood_new-log_likelihood_ref)
                                *gamma_distribution.pdf(d1_new)/gamma_distribution.pdf(d1_old)
                                *gamma_distribution.pdf(d2_new)/gamma_distribution.pdf(d2_old)
                                *jumping_distribution->p(xi_old, xi_new);
                }
                x   = gsl_ran_flat(rng, 0.0, 1.0);
                print_debug("rho: %f\n", rho);
                print_debug("rho: %f\n", std::min(1.0, rho));
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        accept();
                        print_debug("accepted d1 (xi) : %f\n", d1_new);
                        print_debug("accepted d2 (xi) : %f\n", d2_new);
                        update_acceptance(id_right, true);
                        return log_likelihood_new;
                }
                else {
                        // sample rejected
                        reject(id_left, id_right, d1_old, d2_old);
                        print_debug("rejected d1 (xi) : %f\n", d1_new);
                        print_debug("rejected d2 (xi) : %f\n", d2_new);
                        update_acceptance(id_right, false);
                        return log_likelihood_ref;
                }
        }
        double sample(pt_node_t::id_t id, double log_likelihood_ref) {
                log_likelihood_ref = sample_eta(id, log_likelihood_ref);
                log_likelihood_ref = sample_xi (id, log_likelihood_ref);
                print_debug("--------------------------------------------------------------------------------\n");
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
                for (pt_node_t::id_t id = 0; id < nodes; id++) {
                        // if this is not a leaf then sample it
                        if (!trees[0]->node_map[id]->leaf()) {
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

        // for each position of the alignment keep a tree that
        // represents this position and two vectors of polynomials for
        // the likelihood
        typedef std::vector<pt_root_t*> trees_t;
        typedef std::vector<pt_cached_polynomial_t<CODE_TYPE, ALPHABET_SIZE> > polynomials_t;

        pt_node_t::id_t nodes;
        trees_t trees;
        polynomials_t polynomials;
        polynomials_t polynomials_tmp;
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
