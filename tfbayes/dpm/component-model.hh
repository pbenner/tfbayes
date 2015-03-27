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
                : m_cluster_assignments(NULL)
                { }
        component_model_t(const model_id_t& model_id)
                : m_model_id           (model_id)
                , m_cluster_assignments(NULL)
                { }
        component_model_t(const model_id_t& model_id, const data_i<cluster_tag_t>& cluster_assignments)
                : m_model_id           (model_id)
                , m_cluster_assignments(&cluster_assignments)
                { }
        virtual ~component_model_t() { };

        virtual component_model_t* clone() const = 0;

        friend void swap(component_model_t& first, component_model_t& second) {
                using std::swap;
                swap(first.m_model_id,            second.m_model_id);
                swap(first.m_cluster_assignments, second.m_cluster_assignments);
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
        virtual void update(const std::string& msg_prefix = "") { }
        virtual const model_id_t& id() const { return m_model_id; }
        virtual       model_id_t& id()       { return m_model_id; }

        virtual const data_i<cluster_tag_t>& cluster_assignments() const {
                return *m_cluster_assignments;
        }
        virtual void set_cluster_assignments(const data_i<cluster_tag_t>& cluster_assignments) {
                m_cluster_assignments = &cluster_assignments;
        }

protected:
        model_id_t m_model_id;
        const data_i<cluster_tag_t>* m_cluster_assignments;
};

#include <tfbayes/dpm/component-model_foreground.hh>
#include <tfbayes/dpm/component-model_background.hh>
#include <tfbayes/dpm/component-model_gaussian.hh>

#endif /* __TFBAYES_DPM_COMPONENT_MODEL_HH__ */
