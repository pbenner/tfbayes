/* Copyright (C) 2013 Philipp Benner
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

#ifndef __TFBAYES_FG_MESSAGES_HH__
#define __TFBAYES_FG_MESSAGES_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cmath>

#include <tfbayes/fg/distribution.hh>

// a message from a factor node to a variable node
typedef exponential_family_i p_message_t;
// a message from a variable node to a factor node
typedef statistics_t q_message_t;

inline
bool operator==(const q_message_t& rhs, const q_message_t& lhs) {
        if (rhs.size() != lhs.size()) {
                return false;
        }
        for (size_t i = 0; i < rhs.size(); i++) {
                if (std::abs(rhs[i] - lhs[i]) > 1.0e-10) {
                        return false;
                }
        }
        return true;
}

inline
bool operator!=(const q_message_t& rhs, const q_message_t& lhs) {
        return !operator==(rhs, lhs);
}

#endif /* __TFBAYES_FG_MESSAGES_HH__ */
