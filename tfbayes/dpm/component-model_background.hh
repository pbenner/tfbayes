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

#ifndef __TFBAYES_DPM_COMPONENT_MODEL_BACKGROUND_HH__
#define __TFBAYES_DPM_COMPONENT_MODEL_BACKGROUND_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

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
        const sequence_data_t<data_tfbs_t::code_t>& data() const {
                return *_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const independence_background_t& pd);

protected:
        size_t _size;

        cluster_tag_t _bg_cluster_tag;

        sequence_data_t<double> _precomputed_marginal;

        const sequence_data_t<data_tfbs_t::code_t>* _data;
};

// Independence Background Model
////////////////////////////////////////////////////////////////////////////////

class independence_mixture_background_t : public component_model_t {
public:
         independence_mixture_background_t(
                 const std::matrix<double>& alpha,
                 const std::vector<double>& weights,
                 const sequence_data_t<data_tfbs_t::code_t>& data,
                 const sequence_data_t<cluster_tag_t>& cluster_assignments,
                 boost::optional<const alignment_set_t<>&> alignment_set =
                 boost::optional<const alignment_set_t<>&>());
         independence_mixture_background_t(const independence_mixture_background_t& distribution);
        ~independence_mixture_background_t();

        independence_mixture_background_t* clone() const;

        friend void swap(independence_mixture_background_t& first, independence_mixture_background_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first._size,                 second._size);
                swap(first._bg_cluster_tag,       second._bg_cluster_tag);
                swap(first._precomputed_marginal, second._precomputed_marginal);
                swap(first._data,                 second._data);
        }

        independence_mixture_background_t& operator=(const component_model_t& component_model);

        // datatypes
        typedef data_tfbs_t::code_t counts_t;

        void precompute_marginal(
                const std::vector<counts_t>& alpha,
                const std::vector<double  >& weights);
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
        const sequence_data_t<data_tfbs_t::code_t>& data() const {
                return *_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const independence_mixture_background_t& pd);

protected:
        size_t _size;

        cluster_tag_t _bg_cluster_tag;

        sequence_data_t<double> _precomputed_marginal;

        const sequence_data_t<data_tfbs_t::code_t>* _data;
};

// Entropy Background Model
////////////////////////////////////////////////////////////////////////////////

class entropy_background_t : public component_model_t {
public:
         entropy_background_t(
                 const std::vector<double>& parameters,
                 const sequence_data_t<data_tfbs_t::code_t>& data,
                 const sequence_data_t<cluster_tag_t>& cluster_assignments,
                 thread_pool_t& thread_pool,
                 const std::string& cachefile = "",
                 boost::optional<const alignment_set_t<>&> alignment_set =
                 boost::optional<const alignment_set_t<>&>(),
                 size_t verbose = 0);
         entropy_background_t(const entropy_background_t& distribution);
        ~entropy_background_t();

        entropy_background_t* clone() const;

        friend void swap(entropy_background_t& first, entropy_background_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first.m_size,                 second.m_size);
                swap(first.m_bg_cluster_tag,       second.m_bg_cluster_tag);
                swap(first.m_precomputed_marginal, second.m_precomputed_marginal);
                swap(first.m_log_likelihood,       second.m_log_likelihood);
                swap(first.m_data,                 second.m_data);
        }

        entropy_background_t& operator=(const component_model_t& component_model);

        // datatypes
        typedef data_tfbs_t::code_t counts_t;

        bool load_marginal(
                const std::vector<double>& parameters,
                const std::string& cachefile);
        bool save_marginal(
                const std::vector<double>& parameters,
                const std::string& cachefile);
        void precompute_marginal(
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
        const sequence_data_t<data_tfbs_t::code_t>& data() const {
                return *m_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const entropy_background_t& pd);

protected:
        size_t m_size;

        cluster_tag_t m_bg_cluster_tag;

        sequence_data_t<double> m_precomputed_marginal;

        double m_log_likelihood;

        const sequence_data_t<data_tfbs_t::code_t>* m_data;

        size_t m_verbose;
};

// Default Background Model
////////////////////////////////////////////////////////////////////////////////

class default_background_t : public component_model_t {
public:
         default_background_t(
                 const std::matrix<double>& alpha,
                 const std::vector<double>& parameters,
                 const std::vector<double>& weights,
                 const sequence_data_t<data_tfbs_t::code_t>& data,
                 const sequence_data_t<cluster_tag_t>& cluster_assignments,
                 thread_pool_t& thread_pool,
                 const std::string& cachefile = "",
                 boost::optional<const alignment_set_t<>&> alignment_set =
                 boost::optional<const alignment_set_t<>&>(),
                 size_t verbose = 0);
         default_background_t(const default_background_t& distribution);
        ~default_background_t();

        default_background_t* clone() const;

        friend void swap(default_background_t& first, default_background_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first.m_alpha,                 second.m_alpha);
                swap(first.m_weights,               second.m_weights);
                swap(first.m_prior_distribution,    second.m_prior_distribution);
                swap(first.m_size1,                 second.m_size1);
                swap(first.m_size2,                 second.m_size2);
                swap(first.m_bg_cluster_tag,        second.m_bg_cluster_tag);
                swap(first.m_log_likelihood,        second.m_log_likelihood);
                swap(first.m_marginal_probability,  second.m_marginal_probability);
                swap(first.m_component_assignments, second.m_component_assignments);
                swap(first.m_data,                  second.m_data);
                swap(first.m_verbose,               second.m_verbose);
                swap(first.m_g,                     second.m_g);
                swap(first.m_g_prev,                second.m_g_prev);
                swap(first.m_epsilon,               second.m_epsilon);
                swap(first.m_n,                     second.m_n);
        }

        default_background_t& operator=(const component_model_t& component_model);

        // datatypes
        typedef data_tfbs_t::code_t counts_t;

        void update();

        ssize_t max_component(const index_t& index) const;

        bool compute_component_assignments_loop();
        bool compute_component_assignments();
        void compute_marginal();

        void gradient(const index_t& index, size_t i, size_t j, double alpha_sum);
        void gradient(const index_t& index, size_t i, double alpha_sum);
        void gradient();

        bool gradient_ascent();
        bool gradient_ascent_loop(
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
        const sequence_data_t<data_tfbs_t::code_t>& data() const {
                return *m_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const default_background_t& pd);

protected:
        std::vector<counts_t> m_alpha;
        std::vector<double  > m_weights;
        boost::math::gamma_distribution<> m_prior_distribution;

        size_t m_size1;
        size_t m_size2;

        cluster_tag_t m_bg_cluster_tag;

        double m_log_likelihood;
        sequence_data_t<double > m_marginal_probability;
        sequence_data_t<ssize_t> m_component_assignments;
        const sequence_data_t<data_tfbs_t::code_t>* m_data;

        size_t m_verbose;

        /* gradients */
        std::matrix<double> m_g;
        std::matrix<double> m_g_prev;
        std::matrix<double> m_epsilon;
        /* number of positions assigned to each component */
        std::vector<double> m_n;
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
                swap(first.alpha,                   second.alpha);
                swap(first.counts,                  second.counts);
                swap(first.weights,                 second.weights);
                swap(first.m_size1,                 second.m_size1);
                swap(first.m_size2,                 second.m_size2);
                swap(first.m_data,                  second.m_data);
                swap(first.m_component_assignments, second.m_component_assignments);
        }

        mixture_dirichlet_t& operator=(const component_model_t& component_model);

        // datatypes
        typedef data_tfbs_t::code_t counts_t;

        size_t add(const index_t& index);
        size_t add(const range_t& range);
        size_t remove(const index_t& index);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range);
        double log_predictive(const index_t& index);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range);
        double log_likelihood() const;
        std::string print_counts() const;

        const sequence_data_t<data_tfbs_t::code_t>& data() const {
                return *m_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const mixture_dirichlet_t& pd);

protected:
        bool precompute_component_assignments_loop();
        void precompute_component_assignments();
        ssize_t max_component(const index_t& index) const;

        std::vector<counts_t> alpha;
        std::vector<counts_t> counts;
        std::vector<double> weights;

        size_t m_size1;
        size_t m_size2;

        const sequence_data_t<data_tfbs_t::code_t>* m_data;
        sequence_data_t<ssize_t> m_component_assignments;
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


#endif /* __TFBAYES_DPM_COMPONENT_MODEL_BACKGROUND_HH__ */
