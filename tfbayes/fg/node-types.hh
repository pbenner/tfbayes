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

#ifndef FG_NODE_TYPES_HH
#define FG_NODE_TYPES_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <messages.hh>

// It is possible to combine different message passing algorithms in
// one factor graph, e.g. the sum product algorithm for discrete nodes
// with variational message passing for continuous ones. Different
// types of algorithms meet at variable nodes, hence, those nodes need
// to understand messages that originate from different
// algorithms.

class   factor_node_i;
class variable_node_i;

class factor_node_i : public clonable {
public:
        virtual factor_node_i* clone() const = 0;

        // link a variable node to this factor node
        virtual void link(size_t i, variable_node_i* variable_node) = 0;

        // send messages to all connected variable nodes
        virtual void send_messages() = 0;
        // receive a message from a variable node (q message), this
        // method should do nothing but to save the message and notify
        // the factor graph
        virtual void recv_message(const q_message_t& msg) = 0;

protected:
        // send a message to the i'th connected variable
        // node (p messages)
        virtual void send_message(size_t i) = 0;
};

class variable_node_i : public clonable {
public:
        virtual variable_node_i* clone() const = 0;

        // link a factor node to this variable node
        virtual void link(size_t i, factor_node_i* factor_node) = 0;

        // send messages to all connected factor nodes
        virtual void send_messages() = 0;
        // receive a message from a factor node (p message), this
        // method should do nothing but to save the message and notify
        // the factor graph
        virtual void recv_message(const p_message_t& msg) = 0;

protected:
        // send a message to the i'th connected factor
        // node (q messages)
        virtual void send_message(size_t i) = 0;
};

#endif /* FG_NODE_TYPES_HH */
