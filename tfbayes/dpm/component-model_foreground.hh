/* Copyright (C) 2011-2015 Philipp Benner
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

#ifndef __TFBAYES_DPM_COMPONENT_MODEL_FOREGROUND_HH__
#define __TFBAYES_DPM_COMPONENT_MODEL_FOREGROUND_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>
#include <vector>

// Multinomial/Dirichlet Model
////////////////////////////////////////////////////////////////////////////////

class product_dirichlet_t : public component_model_t
{
public:
        template <class S, class T>
        product_dirichlet_t(
                const model_id_t& model_id,
                const S& a_alpha,
                const T& lengths,
                const sequence_data_t<data_tfbs_t::code_t>& data,
                const sequence_data_t<data_tfbs_t::code_t>& complement_data)
                : component_model_t (model_id)
                , m_lengths         (lengths.begin(), lengths.end())
                , m_data            (&data)
                , m_complement_data (&complement_data) {
                // make sure the counts vector has a single column
                assert(a_alpha   .size() == 1);
                assert(a_alpha[0].size() == size2());

                for (size_t i = 0; i < size1(); i++) {
                        m_alpha .push_back(counts_t());
                        m_counts.push_back(counts_t());
                        for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                                m_alpha [i][j]     = a_alpha[0][j];
                                m_counts[i][j]     = a_alpha[0][j];
                                m_alpha_default[j] = a_alpha[0][j];
                        }
                }
                // the lengths should be sorted so that proposals are
                // similar in length
                std::sort(m_lengths.begin(), m_lengths.end());
        }
         product_dirichlet_t(const product_dirichlet_t& distribution);
        ~product_dirichlet_t();

        product_dirichlet_t* clone() const;

        friend void swap(product_dirichlet_t& first, product_dirichlet_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first.m_alpha,           second.m_alpha);
                swap(first.m_counts,          second.m_counts);
                swap(first.m_alpha_default,   second.m_alpha_default);
                swap(first.m_lengths,         second.m_lengths);
                swap(first.m_tmp_counts,      second.m_tmp_counts);
                swap(first.m_data,            second.m_data);
                swap(first.m_complement_data, second.m_complement_data);
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

        const std::vector<size_t>& lengths() const {
                return m_lengths;
        }
        const sequence_data_t<data_tfbs_t::code_t>& data() const {
                return *m_data;
        }
        const sequence_data_t<data_tfbs_t::code_t>& complement_data() const {
                return *m_complement_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const product_dirichlet_t& pd);

protected:
        std::vector<counts_t> m_alpha;
        std::vector<counts_t> m_counts;
        // pseudocounts for resizing the model
        counts_t m_alpha_default;
        // feasible lengths of this model
        std::vector<size_t> m_lengths;

        counts_t m_tmp_counts;

        size_t size1() const { return component_model_t::m_model_id.length; }
        size_t size2() const { return data_tfbs_t::alphabet_size;          }

        const sequence_data_t<data_tfbs_t::code_t>* m_data;
        const sequence_data_t<data_tfbs_t::code_t>* m_complement_data;
};

#endif /* __TFBAYES_DPM_COMPONENT_MODEL_FOREGROUND_HH__ */
