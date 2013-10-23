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

class dpm_subset_t : public clonable {
public:
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
        // typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_set<index_i*, IndexPtrHash, IndexPtrEqual> index_set_t;
        typedef index_set_t::iterator iterator;
        typedef index_set_t::const_iterator const_iterator;

        // constructors
        ////////////////////////////////////////////////////////////////////////
        dpm_subset_t(const dpm_subset_tag_t& dpm_subset_tag)
                : _dpm_subset_tag(dpm_subset_tag)
                {}
        dpm_subset_t(const dpm_subset_t& dpm_subset)
                : _dpm_subset_tag(dpm_subset._dpm_subset_tag) {
                // make sure that every element is cloned
                for (const_iterator it = dpm_subset.begin();
                     it != dpm_subset.end(); it++) {
                        insert(**it);
                }
        }
        // free every index in the set before the set gets destructed
        virtual ~dpm_subset_t() {
                for (index_set_t::const_iterator it = _index_set.begin();
                     it != _index_set.end(); it++) {
                        delete(*it);
                }
        }

        dpm_subset_t* clone() const {
                return new dpm_subset_t(*this);
        }

        friend void swap(dpm_subset_t& left, dpm_subset_t& right) {
                using std::swap;
                swap(left._dpm_subset_tag, right._dpm_subset_tag);
                swap(left._index_set,      right._index_set);
        }

        // operators
        ////////////////////////////////////////////////////////////////////////
        dpm_subset_t& operator=(const dpm_subset_t& dpm_subset) {
                dpm_subset_t tmp(dpm_subset);
                swap(*this, tmp);
                return *this;
        }

        // methods
        ////////////////////////////////////////////////////////////////////////
        std::pair<iterator, bool> insert(const index_i& index) {
                std::pair<iterator, bool> result;
                // clone index
                index_i* tmp = index.clone();
                // insert address to index
                result = _index_set.insert(tmp);
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
                return _index_set.erase(position);
        }
        size_t erase(const index_i& index) {
                // find index pointer
                iterator it = find(index);
                if (it != end()) {
                        // if index exists, erase it
                        erase(it);
                        // one erased
                        return 1;
                }
                // free memory of cloned index
                return 0;
        }
        iterator find(const index_i& index) {
                index_i* i = index.clone();
                // find index pointer
                iterator it = _index_set.find(i);
                // free index
                delete(i);
                return it;
        }
        // each subset represents a cluster with a specific model
        const dpm_subset_tag_t dpm_subset_tag() const {
                return _dpm_subset_tag;
        }
              iterator begin()       { return _index_set.begin(); };
        const_iterator begin() const { return _index_set.begin(); };
              iterator end()         { return _index_set.end();   };
        const_iterator end()   const { return _index_set.end();   };

        size_t size() const {
                return _index_set.size();
        }
        bool operator==(const dpm_subset_t& rhs) const {
                return _dpm_subset_tag == rhs._dpm_subset_tag &&
                        _index_set == rhs._index_set;
        }
        bool operator!=(const dpm_subset_t& rhs) const {
                return !operator==(rhs);
        }
protected:
        dpm_subset_tag_t _dpm_subset_tag;
        index_set_t _index_set;
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
