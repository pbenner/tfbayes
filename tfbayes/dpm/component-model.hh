/* Copyright (C) 2011-2013 Philipp Benner
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

#ifndef __TFBAYES_DPM_COMPONENT_MODEL_HH__
#define __TFBAYES_DPM_COMPONENT_MODEL_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>

#include <boost/math/distributions/gamma.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/utility/clonable.hh>
#include <tfbayes/dpm/data-tfbs.hh>
#include <tfbayes/dpm/datatypes.hh>
#include <tfbayes/dpm/mixture-weights.hh>
#include <tfbayes/dpm/nucleotide-context.hh>
#include <tfbayes/dpm/dpm-tfbs-options.hh>
#include <tfbayes/utility/thread-pool.hh>

// component_model_t interface
////////////////////////////////////////////////////////////////////////////////

class component_model_t : public virtual clonable {
public:
        component_model_t()
                : _cluster_assignments(NULL)
                { }
        component_model_t(const data_i<cluster_tag_t>& cluster_assignments)
                : _cluster_assignments(&cluster_assignments)
                { }
        virtual ~component_model_t() { };

        virtual component_model_t* clone() const = 0;

        friend void swap(component_model_t& first, component_model_t& second) {
                using std::swap;
                swap(first._cluster_assignments, second._cluster_assignments);
        }

        virtual component_model_t& operator=(const component_model_t& component_model) = 0;

        // purely virtual functions
        virtual size_t add(const range_t& range) = 0;
        virtual size_t remove(const range_t& range) = 0;
        virtual size_t count(const range_t& range) = 0;
        virtual double predictive(const range_t& range) = 0;
        virtual double predictive(const std::vector<range_t>& range_set) = 0;
        virtual double log_predictive(const range_t& range) = 0;
        virtual double log_predictive(const std::vector<range_t>& range_set) = 0;
        virtual double log_likelihood() const = 0;
        virtual std::string print_counts() const { return std::string(); }
        virtual void update() { }

        virtual const data_i<cluster_tag_t>& cluster_assignments() const {
                return *_cluster_assignments;
        }
        virtual void set_cluster_assignments(const data_i<cluster_tag_t>& cluster_assignments) {
                _cluster_assignments = &cluster_assignments;
        }

protected:
        const data_i<cluster_tag_t>* _cluster_assignments;
};

// Independence Background Model
////////////////////////////////////////////////////////////////////////////////

class independence_background_t : public component_model_t {
public:
         independence_background_t(
                 const std::vector<double>& alpha,
                 const std::vector<double>& parameters,
                 const sequence_data_t<data_tfbs_t::code_t>& data,
                 const sequence_data_t<cluster_tag_t>& cluster_assignments,
                 thread_pool_t& thread_pool,
                 const std::string& cachefile = "",
                 boost::optional<const alignment_set_t<>&> alignment_set =
                 boost::optional<const alignment_set_t<>&>());
         independence_background_t(const independence_background_t& distribution);
        ~independence_background_t();

        independence_background_t* clone() const;

        friend void swap(independence_background_t& first, independence_background_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first._size,                 second._size);
                swap(first._bg_cluster_tag,       second._bg_cluster_tag);
                swap(first._precomputed_marginal, second._precomputed_marginal);
                swap(first._data,                 second._data);
        }

        independence_background_t& operator=(const component_model_t& component_model);

        // datatypes
        typedef data_tfbs_t::code_t counts_t;

        bool load_marginal_gamma(
                const counts_t& alpha,
                const std::vector<double>& parameters,
                const std::string& cachefile);
        bool save_marginal_gamma(
                const counts_t& alpha,
                const std::vector<double>& parameters,
                const std::string& cachefile);
        void precompute_marginal(
                const counts_t& alpha);
        void precompute_marginal_gamma(
                const counts_t& alpha,
                const std::vector<double>& parameters,
                thread_pool_t& thread_pool);
        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range_set);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range_set);
        double log_likelihood() const;
        std::string print_counts() const;
        void set_bg_cluster_tag(cluster_tag_t cluster_tag);

        const sequence_data_t<cluster_tag_t>& cluster_assignments() const {
                return static_cast<const sequence_data_t<cluster_tag_t>&>(component_model_t::cluster_assignments());
        }
        const sequence_data_t<data_tfbs_t::code_t>& data() {
                return *_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const independence_background_t& pd);

protected:
        size_t _size;

        cluster_tag_t _bg_cluster_tag;

        sequence_data_t<double> _precomputed_marginal;

        const sequence_data_t<data_tfbs_t::code_t>* _data;
};

// Default Background Model
////////////////////////////////////////////////////////////////////////////////

class default_background_t : public component_model_t {
public:
         default_background_t(
                 const std::vector<double>& parameters,
                 const sequence_data_t<data_tfbs_t::code_t>& data,
                 const sequence_data_t<cluster_tag_t>& cluster_assignments,
                 thread_pool_t& thread_pool,
                 const std::string& cachefile = "",
                 boost::optional<const alignment_set_t<>&> alignment_set =
                 boost::optional<const alignment_set_t<>&>());
         default_background_t(const default_background_t& distribution);
        ~default_background_t();

        default_background_t* clone() const;

        friend void swap(default_background_t& first, default_background_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first.alpha,                 second.alpha);
                swap(first.prior_distribution,    second.prior_distribution);
                swap(first._size,                 second._size);
                swap(first._bg_cluster_tag,       second._bg_cluster_tag);
                swap(first._precomputed_marginal, second._precomputed_marginal);
                swap(first._data,                 second._data);
        }

        default_background_t& operator=(const component_model_t& component_model);

        // datatypes
        typedef data_tfbs_t::code_t counts_t;

        void update();
        void precompute_marginal();

        double gradient(const seq_index_t& index, size_t k, double alpha_sum);
        void   gradient(const seq_index_t& index, double alpha_sum, std::vector<double>& result);
        void   gradient(std::vector<double>& result);

        void gradient_ascent();
        double gradient_ascent(
                std::vector<double>& g,
                std::vector<double>& g_prev,
                std::vector<double>& epsilon,
                double eta = 0.1,
                double min_alpha = 1.0e-20);

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range_set);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range_set);
        double log_likelihood() const;
        std::string print_counts() const;
        void set_bg_cluster_tag(cluster_tag_t cluster_tag);

        const sequence_data_t<cluster_tag_t>& cluster_assignments() const {
                return static_cast<const sequence_data_t<cluster_tag_t>&>(component_model_t::cluster_assignments());
        }
        const sequence_data_t<data_tfbs_t::code_t>& data() {
                return *_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const default_background_t& pd);

protected:
        counts_t alpha;
        boost::math::gamma_distribution<> prior_distribution;

        size_t _size;

        cluster_tag_t _bg_cluster_tag;

        sequence_data_t<double> _precomputed_marginal;
        const sequence_data_t<data_tfbs_t::code_t>* _data;
};

// Multinomial/Dirichlet Model
////////////////////////////////////////////////////////////////////////////////

class product_dirichlet_t : public component_model_t {
public:
         product_dirichlet_t(const std::matrix<double>& alpha,
                             const sequence_data_t<data_tfbs_t::code_t>& data,
                             const sequence_data_t<data_tfbs_t::code_t>& complement_data,
                             const bool complement = true);
         product_dirichlet_t(const product_dirichlet_t& distribution);
        ~product_dirichlet_t();

        product_dirichlet_t* clone() const;

        friend void swap(product_dirichlet_t& first, product_dirichlet_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first.alpha,            second.alpha);
                swap(first.counts,           second.counts);
                swap(first.tmp_counts,       second.tmp_counts);
                swap(first._size1,           second._size1);
                swap(first._size2,           second._size2);
                swap(first._data,            second._data);
                swap(first._complement_data, second._complement_data);
                swap(first._complement,      second._complement);
        }

        product_dirichlet_t& operator=(const component_model_t& component_model);

        // datatypes
        typedef data_tfbs_t::code_t counts_t;

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range);
        double log_likelihood() const;
        std::string print_counts() const;

        const sequence_data_t<data_tfbs_t::code_t>& data() const {
                return *_data;
        }
        const sequence_data_t<data_tfbs_t::code_t>& complement_data() const {
                return *_complement_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const product_dirichlet_t& pd);

protected:
        std::vector<counts_t> alpha;
        std::vector<counts_t> counts;

        counts_t tmp_counts;

        size_t _size1;
        size_t _size2;

        const sequence_data_t<data_tfbs_t::code_t>* _data;
        const sequence_data_t<data_tfbs_t::code_t>* _complement_data;
        bool _complement;
};

// Multinomial/Dirichlet Mixture Model
////////////////////////////////////////////////////////////////////////////////

class mixture_dirichlet_t : public component_model_t {
public:
         mixture_dirichlet_t(const std::matrix<double>& alpha,
                             const std::vector<double>& weights,
                             const sequence_data_t<data_tfbs_t::code_t>& data);
         mixture_dirichlet_t(const mixture_dirichlet_t& distribution);
        ~mixture_dirichlet_t();

        mixture_dirichlet_t* clone() const;

        friend void swap(mixture_dirichlet_t& first, mixture_dirichlet_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first.alpha,            second.alpha);
                swap(first.counts,           second.counts);
                swap(first.weights,          second.weights);
                swap(first._size1,           second._size1);
                swap(first._size2,           second._size2);
                swap(first._data,            second._data);
                swap(first._component_assignments, second._component_assignments);
        }

        mixture_dirichlet_t& operator=(const component_model_t& component_model);

        // datatypes
        typedef data_tfbs_t::code_t counts_t;

        size_t add(const index_i& index);
        size_t add(const range_t& range);
        size_t remove(const index_i& index);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range);
        double log_predictive(const index_i& index);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range);
        double log_likelihood() const;
        std::string print_counts() const;

        const sequence_data_t<data_tfbs_t::code_t>& data() const {
                return *_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const mixture_dirichlet_t& pd);

protected:
        bool precompute_component_assignments_loop();
        void precompute_component_assignments();
        ssize_t max_component(const seq_index_t& index) const;

        std::vector<counts_t> alpha;
        std::vector<counts_t> counts;
        std::vector<double> weights;

        size_t _size1;
        size_t _size2;

        const sequence_data_t<data_tfbs_t::code_t>* _data;
        sequence_data_t<ssize_t> _component_assignments;
};

// Markov Chain Mixture
////////////////////////////////////////////////////////////////////////////////

class markov_chain_mixture_t : public component_model_t {
public:
         markov_chain_mixture_t(size_t alphabet_size,
                                const tfbs_options_t& options,
                                const sequence_data_t<data_tfbs_t::code_t>& data,
                                const sequence_data_t<cluster_tag_t>& cluster_assignments,
                                cluster_tag_t cluster_tag);
         markov_chain_mixture_t(const markov_chain_mixture_t& distribution);
        ~markov_chain_mixture_t();

        markov_chain_mixture_t* clone() const;

        markov_chain_mixture_t& operator=(const component_model_t& component_model);

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range_set);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range_set);
        double log_likelihood() const;

        const sequence_data_t<cluster_tag_t>& cluster_assignments() const {
                return static_cast<const sequence_data_t<cluster_tag_t>&>(component_model_t::cluster_assignments());
        }

        friend std::ostream& operator<< (std::ostream& o, const markov_chain_mixture_t& pd);

protected:
        size_t  _length;
        double* _counts;
        double* _alpha;
        /* sum of pseudocounts and real counts for each character
         * followed by a specific context */
        double* _counts_sum;
        /* temporary copy of counts and counts_sum to compute
         * likelihoods */
        double* _counts_tmp;
        int * _parents;

        mixture_weights_t* _weights;

        const sequence_data_t<data_tfbs_t::code_t>& _data;
              sequence_data_t<context_t>            _context;

        const cluster_tag_t _cluster_tag;
        const size_t        _max_context;
        const size_t        _alphabet_size;

        // internal methods
        size_t max_from_context(const range_t& range) const;
        size_t max_to_context(const range_t& range) const;
        double log_likelihood(size_t pos) const;
        void substract_counts(size_t pos) const;
};

// Bivariate Gaussian
////////////////////////////////////////////////////////////////////////////////

class bivariate_normal_t : public component_model_t {
public:
         bivariate_normal_t();
         bivariate_normal_t(const std::matrix<double>& Sigma,
                            const std::matrix<double>& Sigma_0,
                            const std::vector<double>& mu_0,
                            const data_t<std::vector<double> >& data);
         bivariate_normal_t(const bivariate_normal_t& bn);
        ~bivariate_normal_t();

        bivariate_normal_t* clone() const;

        friend void swap(bivariate_normal_t& first, bivariate_normal_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first._Sigma_0_inv, second._Sigma_0_inv);
                swap(first._mu_0,        second._mu_0);
                swap(first._Sigma,       second._Sigma);
                swap(first._Sigma_inv,   second._Sigma_inv);
                swap(first._mu,          second._mu);
                swap(first._N,           second._N);
                swap(first._Sigma_N,     second._Sigma_N);
                swap(first._Sigma_N_inv, second._Sigma_N_inv);
                swap(first._mu_N,        second._mu_N);
                swap(first._dimension,   second._dimension);
                swap(first._inv_perm,    second._inv_perm);
                swap(first._inv_tmp,     second._inv_tmp);
                swap(first._tmp1,        second._tmp1);
                swap(first._tmp2,        second._tmp2);
                swap(first._data,        second._data);
        }

        bivariate_normal_t& operator=(const component_model_t& component_model);

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range_set);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range_set);
        double log_likelihood() const;
        std::vector<double> mean() const;

        const data_t<std::vector<double> >& data() const {
                return *_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const bivariate_normal_t& pd);

protected:
        // prior
        gsl_matrix* _Sigma_0_inv;
        gsl_vector* _mu_0;

        // likelihood
        gsl_matrix* _Sigma;
        gsl_matrix* _Sigma_inv;
        gsl_vector* _mu;
        double _N;

        // posterior
        gsl_matrix* _Sigma_N;
        gsl_matrix* _Sigma_N_inv;
        gsl_vector* _mu_N;

        // other
        size_t _dimension;
        gsl_permutation* _inv_perm;
        gsl_matrix* _inv_tmp;
        gsl_vector* _tmp1;
        gsl_vector* _tmp2;

        void inverse(gsl_matrix* dst, const gsl_matrix* src);
        void update();

        const data_t<std::vector<double> >* _data;
};

#endif /* __TFBAYES_DPM_COMPONENT_MODEL_HH__ */
