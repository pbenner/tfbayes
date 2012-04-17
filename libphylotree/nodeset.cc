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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <nodeset.hh>

using namespace std;

ostream& operator<< (ostream& o, const nodeset_t& nodeset) {
        for (nodeset_t::const_iterator it = nodeset.begin(); it != nodeset.end(); it++) {
                if (it != nodeset.begin()) {
                        o << ",";
                }
                o << (*it)->name;
        }

        return o;
}

size_t hash_value(const nodeset_t& set) {
        size_t seed = 0;
        for (nodeset_t::const_iterator it = set.begin(); it != set.end(); it++) {
                seed += (size_t)*it;
        }
        return seed;
}
