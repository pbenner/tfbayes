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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/dpm/index.hh>

#include <boost/format.hpp>

std::ostream&
operator<< (std::ostream& o, const index_t& index) {
        if (index.m_x[1] == -1) {
                o << boost::format("(%d)")
                        % index.m_x[0];
        }
        else {
                o << boost::format("(%d, %d)")
                        % index.m_x[0] % index.m_x[1];
        }
        return o;
}

std::ostream&
operator<< (std::ostream& o, const range_t& range) {
        if (typeid(range.index()) == typeid(index_t)) {
                o << range.index()
                  << ":" << range.length();
        }
        else {
                o << range.index()
                  << ":" << range.length();
        }
        if (range.reverse()) {
                o << "!";
        }
        return o;
}
