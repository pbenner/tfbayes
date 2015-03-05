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

#ifndef __TFBAYES_DPM_DPM_PARTITION_HH__
#define __TFBAYES_DPM_DPM_PARTITION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <utility> // std::pair

#include <boost/unordered_set.hpp>

#include <tfbayes/dpm/index.hh>
#include <tfbayes/dpm/datatypes.hh>

class dpm_subset_t : public virtual clonable, public boost::unordered_set<range_t> {
public:
        // typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_set<range_t> base_t;
        typedef base_t::iterator iterator;
        typedef base_t::const_iterator const_iterator;

        // constructors
        ////////////////////////////////////////////////////////////////////////
        dpm_subset_t(const model_id_t& model_id)
                : base_t()
                , _model_id(model_id)
                {}
        dpm_subset_t(const dpm_subset_t& dpm_subset)
                : base_t(dpm_subset),
                  _model_id(dpm_subset._model_id) {
        }
        virtual ~dpm_subset_t() {
        }

        dpm_subset_t* clone() const {
                return new dpm_subset_t(*this);
        }

        friend void swap(dpm_subset_t& left, dpm_subset_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
                swap(left._model_id, right._model_id);
        }

        // operators
        ////////////////////////////////////////////////////////////////////////
        dpm_subset_t& operator=(const dpm_subset_t& dpm_subset) {
                using std::swap;
                dpm_subset_t tmp(dpm_subset);
                swap(*this, tmp);
                return *this;
        }
        // methods
        ////////////////////////////////////////////////////////////////////////
        // each subset represents a cluster with a specific model
        const model_id_t model_id() const {
                return _model_id;
        }
protected:
        model_id_t _model_id;
};

class dpm_partition_t : public std::vector<dpm_subset_t> {
public:
        dpm_partition_t()
                : std::vector<dpm_subset_t>()
                {}
        template <typename InputIterator>
        dpm_partition_t (InputIterator first, InputIterator last)
                : std::vector<dpm_subset_t>(first, last) { }

        void add_component(const model_id_t& model_id) {
                this->push_back(dpm_subset_t(model_id));
        }
};

typedef std::vector<dpm_partition_t> dpm_partition_list_t;

std::ostream& operator<< (std::ostream& o, const dpm_subset_t& dpm_subset);
std::ostream& operator<< (std::ostream& o, const dpm_partition_t& partition);
std::ostream& operator<< (std::ostream& o, const dpm_partition_list_t& partitions);

#endif /* __TFBAYES_DPM_DPM_PARTITION_HH__ */
