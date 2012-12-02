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

#include <index.hh>

class index_set_t : public boost::unordered_set<index_i*> {
public:
        // typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_set<index_i*> set_t;

        // constructors
        ////////////////////////////////////////////////////////////////////////
        index_set_t()
                : boost::unordered_set<index_i*>()
                {}
        index_set_t(const index_set_t& index_set)
                : boost::unordered_set<index_i*>() {
                operator=(index_set);
        }
        // free every index in the set before the set gets destructed
        virtual ~index_set_t() {
                for (set_t::const_iterator it = this->begin();
                     it != this->end(); it++) {
                        index_i* index = *it;
                        delete(index);
                }
        }
        // make sure that every element is cloned
        index_set_t& operator=(const index_set_t& index_set) {
                this->clear();
                for (set_t::const_iterator it = index_set.begin();
                     it != index_set.end(); it++) {
                        insert(**it);
                }

                return *this;
        }

        // methods
        ////////////////////////////////////////////////////////////////////////
        void insert(const index_i& index) {

                set_t::insert(index.clone());
        }
};

class dpm_partition_t : public std::vector<index_set_t> {
public:
        dpm_partition_t()
                : std::vector<index_set_t>()
                {}

        void add_component() {
                this->push_back(index_set_t());
        }
        
};

#endif /* DPM_PARTITION_HH */
