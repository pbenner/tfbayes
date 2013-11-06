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

#ifndef __TFBAYES_FG_NODE_TYPES_HH__
#define __TFBAYES_FG_NODE_TYPES_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>    // std::fill
#include <cassert>

#include <boost/function.hpp>
#include <boost/bind.hpp>
// std::mutex does not work with Mac OS X, so instead
// we use boost::mutex
#include <boost/thread.hpp>

#include <mailbox.hh>
#include <messages.hh>
#include <observable.hh>

// It is possible to combine different message passing algorithms in
// one factor graph, e.g. the sum product algorithm for discrete nodes
// with variational message passing for continuous ones. Different
// types of algorithms meet at variable nodes, hence, those nodes need
// to understand messages that originate from different
// algorithms and multiply them to get a new q-message. The common
// language, i.e. the typ of q-messages, is exponential families.

// interfaces
////////////////////////////////////////////////////////////////////////////////

// a general node in a factor graph must be clonable, observable,
// and be able to send messages to neighboring nodes
class factor_graph_node_i : public virtual clonable, public virtual observable_i {
public:
        virtual ~factor_graph_node_i() { }

        // send messages to all connected variable nodes
        virtual void send_messages() = 0;
};

class   factor_node_i;
class variable_node_i;

class factor_node_i : public virtual factor_graph_node_i {
public:
        virtual ~factor_node_i() { }

        virtual factor_node_i* clone() const = 0;

        // link a variable node to this factor node
        virtual void link(size_t i, variable_node_i& variable_node) = 0;

        // neighboring variable nodes
        virtual const std::vector<variable_node_i*>& neighbors() const = 0;

protected:
        // prepare a message to the i'th connected variable
        // node (p messages)
        virtual const p_message_t& message(size_t i) = 0;
};

class variable_node_i : public virtual factor_graph_node_i {
public:
        virtual ~variable_node_i() { }

        virtual variable_node_i* clone() const = 0;

        // send messages to all connected factor nodes
        virtual void send_messages() = 0;

        // link a factor node to this variable node
        virtual mailbox_slot_t<p_message_t>& link(mailbox_slot_t<q_message_t>& slot) = 0;
};

// basic implementations
////////////////////////////////////////////////////////////////////////////////

template <size_t D>
class factor_node_t : public virtual factor_node_i, public observable_t {
public:
        factor_node_t() :
                _inbox    (D),
                outbox    (D),
                _neighbors(D, NULL) {
        }
        factor_node_t(const factor_node_t& factor_node) :
        // do not copy the mailer and mailbox, since they should be
        // populated manually to create a new network
                _inbox    (D),
                outbox    (D),
                _neighbors(D, NULL) {
        }

        virtual factor_node_t* clone() const = 0;

        friend void swap(factor_node_t& left, factor_node_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left._inbox,     right._inbox);
                swap(left.outbox,     right.outbox);
                swap(left._neighbors, right._neighbors);
        }

        factor_node_t& operator=(const factor_node_t& node) {
                using std::swap;
                factor_node_t tmp(node);
                swap(*this, tmp);
                return *this;
        }

        virtual void send_messages() {
                std::cout << "factor node " << this << " is sending messages" << std::endl;
                lock_inbox();
                for (size_t i = 0; i < D; i++) {
                        if (outbox[i]) {
                                outbox[i]->lock();
                                outbox[i]->replace(message(i));
                                outbox[i]->notify();
                                outbox[i]->unlock();
                        }
                }
                unlock_inbox();
        }
        virtual void link(size_t i, variable_node_i& variable_node) {
                assert(i < D);
                // pointer to the method that notifies the factor
                // graph about an update
                void (observable_t::*tmp) () const = &factor_node_t::notify;
                // prepare a new inbox slot
                _inbox[i] = *new mailbox_slot_t<q_message_t>(boost::bind(tmp, this));
                // receive and send mailer
                outbox[i] = variable_node.link(*_inbox[i]);
                // initialize outbox
                outbox[i]->replace(initial_message(i));
                // save neighbor
                _neighbors[i] = &variable_node;
        }
        virtual const std::vector<variable_node_i*>& neighbors() const {
                return _neighbors;
        }

protected:
        // initialize outbox i
        virtual const p_message_t& initial_message(size_t i) const = 0;
        // mailboxes
        _inbox_t<q_message_t> _inbox;
        outbox_t<p_message_t> outbox;
        // keep track of neighboring nodes for cloning whole networks
        std::vector<variable_node_i*> _neighbors;
private:
        // lock every slot in the mailbox, so that we can prepare a
        // new message without receiving new mail while doing so
        void lock_inbox() const {
                for (size_t i = 0; i < D; i++) {
                        if (_inbox[i]) _inbox[i]->lock();
                }
        }
        void unlock_inbox() const {
                for (size_t i = 0; i < D; i++) {
                        if (_inbox[i]) _inbox[i]->unlock();
                }
        }
};

template <typename T>
class variable_node_t : public virtual variable_node_i, public observable_t {
public:
        variable_node_t() :
                _inbox      (),
                outbox      (),
                old_message (),
                new_message () {
        }
        variable_node_t(const variable_node_t& variable_node) :
        // do not copy the inbox and outbox, since they should be
        // populated manually to create a new network
                _inbox      (),
                outbox      (),
                old_message (variable_node.old_message),
                new_message (variable_node.new_message) {
        }
        virtual variable_node_t* clone() const {
                return new variable_node_t(*this);
        }

        friend void swap(variable_node_t& left, variable_node_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left._inbox,      right._inbox);
                swap(left.outbox,      right.outbox);
                swap(left.old_message, right.old_message);
                swap(left.new_message, right.new_message);
        }

        variable_node_t& operator=(const variable_node_t& node) {
                using std::swap;
                variable_node_t tmp(node);
                swap(*this, tmp);
                return *this;
        }

        virtual void send_messages() {
                std::cout << "variable node " << this << " is sending messages" << std::endl;
                // a temporary message
                T msg;
                // loop over all slots of the mailbox
                for (size_t i = 0; i < _inbox.size(); i++) {
                        // lock this slot
                        _inbox[i]->lock();
                        // and get the message
                        std::cout << "-> getting message from slot " << i << std::endl;
                        msg *= static_cast<const T&>(_inbox[i]->receive());
                        // release lock
                        _inbox[i]->unlock();
                }
                // check if this message was sent before
                if (msg == old_message) return;
                // lock all connected nodes
                for (size_t i = 0; i < outbox.size(); i++) {
                        std::cout << "-> sending message to neighbor " << i << std::endl;
                        outbox[i]->lock();
                }
                // send message
                new_message = msg;
                // unlock and notify
                for (size_t i = 0; i < outbox.size(); i++) {
                        outbox[i]->notify();
                        outbox[i]->unlock();
                }
                // save message
                old_message = new_message;
        }
        virtual mailbox_slot_t<p_message_t>& link(mailbox_slot_t<q_message_t>& slot) {
                void (observable_t::*tmp) () const = &variable_node_t::notify;
                // save slot to the outbox
                outbox.push_back(slot);
                // put the current message into the box
                slot.replace(new_message);
                // and prepare a new inbox for this node
                _inbox.push_back(*new mailbox_slot_t<p_message_t>(boost::bind(tmp, this)));
                return *_inbox[_inbox.size()-1];
        }

protected:
        // mailboxes
        _inbox_t<p_message_t> _inbox;
        outbox_t<q_message_t> outbox;

        T old_message;
        T new_message;
};

#endif /* __TFBAYES_FG_NODE_TYPES_HH__ */
