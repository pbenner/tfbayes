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

// Multinomial/Dirichlet Model
////////////////////////////////////////////////////////////////////////////////

class product_dirichlet_t : public component_model_t {
public:
         product_dirichlet_t(const model_id_t& model_id,
                             const std::matrix<double>& alpha,
                             const sequence_data_t<data_tfbs_t::code_t>& data,
                             const sequence_data_t<data_tfbs_t::code_t>& complement_data);
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
        size_t length() const;

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
};

#endif /* __TFBAYES_DPM_COMPONENT_MODEL_FOREGROUND_HH__ */
