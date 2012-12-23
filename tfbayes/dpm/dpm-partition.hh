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

#ifndef DPM_PARTITION_HH
#define DPM_PARTITION_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/unordered_set.hpp>

#include <tfbayes/dpm/index.hh>
#include <tfbayes/dpm/datatypes.hh>

typedef ssize_t dpm_subset_tag_t;

class dpm_subset_t : public boost::unordered_set<index_i*> {
public:
        // typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_set<index_i*> set_t;

        // constructors
        ////////////////////////////////////////////////////////////////////////
        dpm_subset_t(const dpm_subset_tag_t& dpm_subset_tag)
        : boost::unordered_set<index_i*>(),
          _dpm_subset_tag(dpm_subset_tag)
                {}
        dpm_subset_t(const dpm_subset_t& dpm_subset)
                : boost::unordered_set<index_i*>() {
                operator=(dpm_subset);
        }
        // free every index in the set before the set gets destructed
        virtual ~dpm_subset_t() {
                for (set_t::const_iterator it = this->begin();
                     it != this->end(); it++) {
                        index_i* index = *it;
                        delete(index);
                }
        }
        // make sure that every element is cloned
        dpm_subset_t& operator=(const dpm_subset_t& dpm_subset) {
                this->clear();
                for (set_t::const_iterator it = dpm_subset.begin();
                     it != dpm_subset.end(); it++) {
                        insert(**it);
                }
                _dpm_subset_tag = dpm_subset.dpm_subset_tag();
                return *this;
        }

        // methods
        ////////////////////////////////////////////////////////////////////////
        void insert(const index_i& index) {

                set_t::insert(index.clone());
        }

        // each subset represents a cluster with a specific model
        const dpm_subset_tag_t dpm_subset_tag() const {
                return _dpm_subset_tag;
        }

private:
        dpm_subset_tag_t _dpm_subset_tag;
};

class dpm_partition_t : public std::vector<dpm_subset_t> {
public:
        dpm_partition_t()
                : std::vector<dpm_subset_t>()
                {}

        void add_component(const dpm_subset_tag_t& dpm_subset_tag) {
                this->push_back(dpm_subset_t(dpm_subset_tag));
        }
        
};

#endif /* DPM_PARTITION_HH */
