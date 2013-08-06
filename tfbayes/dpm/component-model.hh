/* Copyright (C) 2011 Philipp Benner
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

#ifndef COMPONENT_MODEL_HH
#define COMPONENT_MODEL_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>

#include <tfbayes/utility/clonable.hh>
#include <tfbayes/dpm/data-tfbs.hh>
#include <tfbayes/dpm/datatypes.hh>
#include <tfbayes/dpm/mixture-weights.hh>
#include <tfbayes/dpm/nucleotide-context.hh>
#include <tfbayes/dpm/dpm-tfbs-options.hh>

// component_model_t interface
////////////////////////////////////////////////////////////////////////////////

class component_model_t : public clonable {
public:
        component_model_t()
                : _cluster_assignments(NULL)
                { }
        component_model_t(const data_t<cluster_tag_t>& cluster_assignments)
                : _cluster_assignments(&cluster_assignments)
                { }

        virtual component_model_t* clone() const = 0;

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

        virtual const data_t<cluster_tag_t>& cluster_assignments() const {
                return *_cluster_assignments;
        }
        virtual void set_cluster_assignments(const data_t<cluster_tag_t>& cluster_assignments) {
                _cluster_assignments = &cluster_assignments;
        }

protected:
        const data_t<cluster_tag_t>* _cluster_assignments;
};

// Independence Background Model
////////////////////////////////////////////////////////////////////////////////

class independence_background_t : public component_model_t {
public:
         independence_background_t(
                 const std::matrix<double>& alpha,
                 const sequence_data_t<data_tfbs_t::code_t>& data,
                 const sequence_data_t<cluster_tag_t>& cluster_assignments);
         independence_background_t(
                 const double k, const double g,
                 const sequence_data_t<data_tfbs_t::code_t>& data,
                 const sequence_data_t<cluster_tag_t>& cluster_assignments);
         independence_background_t(const independence_background_t& distribution);
        ~independence_background_t();

        independence_background_t* clone() const;

        // datatypes
        typedef data_tfbs_t::code_t counts_t;

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range_set);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range_set);
        double log_likelihood() const;
        virtual std::string print_counts() const;
        void set_bg_cluster_tag(cluster_tag_t cluster_tag);

        virtual const sequence_data_t<cluster_tag_t>& cluster_assignments() const {
                return *static_cast<const sequence_data_t<cluster_tag_t>*>(&component_model_t::cluster_assignments());
        }

        friend std::ostream& operator<< (std::ostream& o, const independence_background_t& pd);

protected:
        const sequence_data_t<data_tfbs_t::code_t>& _data;

        const size_t _size;

        cluster_tag_t _bg_cluster_tag;

        sequence_data_t<double> _precomputed_marginal;
};

// Multinomial/Dirichlet Model
////////////////////////////////////////////////////////////////////////////////

class product_dirichlet_t : public component_model_t {
public:
         product_dirichlet_t(const std::matrix<double>& alpha,
                             const sequence_data_t<data_tfbs_t::code_t>& data);
         product_dirichlet_t(const product_dirichlet_t& distribution);
        ~product_dirichlet_t();

        product_dirichlet_t* clone() const;

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
        virtual std::string print_counts() const;

        friend std::ostream& operator<< (std::ostream& o, const product_dirichlet_t& pd);

protected:
        std::vector<counts_t> alpha;
        std::vector<counts_t> counts;

        counts_t tmp_counts;

        const sequence_data_t<data_tfbs_t::code_t>& _data;

        const size_t _size1;
        const size_t _size2;
};

// Void Model
////////////////////////////////////////////////////////////////////////////////

class uniform_background_t : public component_model_t {
public:
         uniform_background_t();
         uniform_background_t(const uniform_background_t& distribution);
        ~uniform_background_t();

        uniform_background_t* clone() const;

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
        virtual std::string print_counts() const;

        friend std::ostream& operator<< (std::ostream& o, const uniform_background_t& pd);

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

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range_set);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range_set);
        double log_likelihood() const;

        virtual const sequence_data_t<cluster_tag_t>& cluster_assignments() const {
                return *static_cast<const sequence_data_t<cluster_tag_t>*>(&component_model_t::cluster_assignments());
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
         bivariate_normal_t(const gsl_matrix* Sigma,
                            const gsl_matrix* Sigma_0,
                            const gsl_vector* mu_0,
                            const data_t<std::vector<double> >& data);
         bivariate_normal_t(const bivariate_normal_t& bn);
        ~bivariate_normal_t();

        bivariate_normal_t* clone() const;

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range_set);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range_set);
        double log_likelihood() const;
        const gsl_vector* mean() const;

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
        const size_t _dimension;
        gsl_permutation* _inv_perm;
        gsl_matrix* _inv_tmp;
        gsl_vector* _tmp1;
        gsl_vector* _tmp2;

        void inverse(gsl_matrix* dst, const gsl_matrix* src);
        void update();

        const data_t<std::vector<double> >& _data;
};

#endif /* COMPONENT_MODEL_HH */
