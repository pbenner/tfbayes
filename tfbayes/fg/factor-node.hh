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

#ifndef FG_FACTOR_NODE_HH
#define FG_FACTOR_NODE_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <distribution.hh>

// Tt is possible to combine different message passing algorithms in
// one factor graph, e.g. the sum product algorithm for discrete nodes
// with variational message passing for continuous ones. Different
// types of algorithms meet at variable nodes, hence, those nodes need
// to understand messages that originate from different
// algorithms. However, variable nodes merely compute the product of
// messages and send them back to adjacent factor nodes.

// a message from a factor node to a variable node
class p_message_i {};
// a message from a variable node to a factor node
class q_message_i {};

// a factor node is an exponential family with the ability to send
// messages to connected nodes
class factor_node_i : public exponential_family_i {
public:
        // send a message to the i'th connected variable
        // node (p messages)
        const p_message_i& send_message(size_t i) const = 0;
        // receive a message from a variable node (q message), this
        // method shoud do nothing but to save the message
        void recv_message(const q_message_i& msg) = 0;
};

class variational_factor_node_i : public factor_node_i, public exponential_family_i {
public:
        // send a message to the i'th connected variable
        // node (p messages)
        const p_message_i& send_message(size_t i) const = 0;
        // receive a message from a variable node (q message), this
        // method shoud do nothing but to save the message
        void recv_message(const q_message_i& msg) = 0;
};

#endif /* FG_FACTOR_NODE_HH */
