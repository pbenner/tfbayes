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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/unordered_set.hpp>

#include <tfbayes/dpm/index.hh>
#include <tfbayes/dpm/datatypes.hh>

typedef std::string dpm_subset_tag_t;

struct IndexPtrHash {
        size_t operator()(index_i* const& index) const {
                return index->hash();
        }
};

struct IndexPtrEqual {
        bool operator()(index_i* const& lhs, index_i* const& rhs) const {
                return *lhs == *rhs;
        }
};

class dpm_subset_t : public boost::unordered_set<index_i*, IndexPtrHash, IndexPtrEqual> {
public:

        // typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_set<index_i*, IndexPtrHash, IndexPtrEqual> base_t;

        // constructors
        ////////////////////////////////////////////////////////////////////////
        dpm_subset_t(const dpm_subset_tag_t& dpm_subset_tag)
                : base_t(),
          _dpm_subset_tag(dpm_subset_tag)
                {}
        dpm_subset_t(const dpm_subset_t& dpm_subset)
                : base_t() {
                operator=(dpm_subset);
        }
        // free every index in the set before the set gets destructed
        virtual ~dpm_subset_t() {
                for (base_t::const_iterator it = this->begin();
                     it != this->end(); it++) {
                        index_i* index = *it;
                        delete(index);
                }
        }

        // operators
        ////////////////////////////////////////////////////////////////////////
        // make sure that every element is cloned
        dpm_subset_t& operator=(const dpm_subset_t& dpm_subset) {
                this->clear();
                for (base_t::const_iterator it = dpm_subset.begin();
                     it != dpm_subset.end(); it++) {
                        insert(**it);
                }
                _dpm_subset_tag = dpm_subset.dpm_subset_tag();
                return *this;
        }
        // methods
        ////////////////////////////////////////////////////////////////////////
        std::pair<iterator, bool> insert(const index_i& index) {
                std::pair<iterator, bool> result;
                // clone index
                index_i* tmp = index.clone();
                // insert address to index
                result = base_t::insert(tmp);
                // free index if insert was not successful
                if (!result.second) {
                        delete(tmp);
                }
                return result;
        }
        iterator erase(const_iterator position) {
                // free index memory
                delete(*position);
                // erase index pointer
                return base_t::erase(position);
        }
        size_type erase(const index_i& index) {
                size_type result = 0;
                // clone the index
                index_i* i = index.clone();
                // find index pointer
                iterator it = find(i);
                if (it != end()) {
                        // if index exists, erase it
                        index_i* tmp = *it;
                        result = base_t::erase(i);
                        // and free memory
                        delete(tmp);
                }
                // free memory of cloned index
                delete(i);
                return result;
        }
        // each subset represents a cluster with a specific model
        const dpm_subset_tag_t dpm_subset_tag() const {
                return _dpm_subset_tag;
        }

protected:
        dpm_subset_tag_t _dpm_subset_tag;
};

class dpm_partition_t : public std::vector<dpm_subset_t> {
public:
        dpm_partition_t()
                : std::vector<dpm_subset_t>()
                {}
        template <typename InputIterator>
        dpm_partition_t (InputIterator first, InputIterator last)
                : std::vector<dpm_subset_t>(first, last) { }

        void add_component(const dpm_subset_tag_t& dpm_subset_tag) {
                this->push_back(dpm_subset_t(dpm_subset_tag));
        }
};

typedef std::vector<dpm_partition_t> dpm_partition_list_t;

std::ostream& operator<< (std::ostream& o, const dpm_subset_t& dpm_subset);
std::ostream& operator<< (std::ostream& o, const dpm_partition_t& partition);
std::ostream& operator<< (std::ostream& o, const dpm_partition_list_t& partitions);

#endif /* DPM_PARTITION_HH */
