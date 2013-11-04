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

#ifndef FG_VARIATIONAL_FACTOR_NODE_HH
#define FG_VARIATIONAL_FACTOR_NODE_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <node-types.hh>

// a factor node is an exponential family with the ability to send
// and receive messages
class normal_node_t : public factor_node_t<3> {
public:
        typedef factor_node_t<3> base_t;

        virtual normal_node_t* clone() const {
                return new normal_node_t(*this);
        }

protected:
        const p_message_t& message(size_t i) {
                return msg1;
        }

        normal_distribution_t distribution;

        normal_distribution_t msg1;
};

#endif /* FG_VARIATIONAL_FACTOR_NODE_HH */
