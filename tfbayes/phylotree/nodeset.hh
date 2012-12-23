/* Copyright (C) 2012 Philipp Benner
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

#ifndef NODESET_HH
#define NODESET_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/unordered_set.hpp>

#include <tfbayes/phylotree/phylotree.hh>

class nodeset_t : public boost::unordered_set<const pt_node_t*> {
public:
        bool empty() const {
                return size() == 0;
        }
        void join(const nodeset_t& nodeset) {
                for (nodeset_t::const_iterator it = nodeset.begin(); it != nodeset.end(); it++) {
                        insert(*it);
                }
        }
};

#include <ostream>

size_t hash_value(const nodeset_t& set);

std::ostream& operator<< (std::ostream& o, const nodeset_t& nodeset);

#endif /* NODESET_HH */
