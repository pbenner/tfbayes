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
#include <typeinfo>
#include <cassert>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/type_traits/is_class.hpp>
#include <boost/utility/enable_if.hpp>
// std::mutex does not work with Mac OS X, so instead
// we use boost::mutex
#include <boost/thread.hpp>

#include <hotnews.hh>
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
        virtual bool link(size_t i, variable_node_i& variable_node) = 0;
        virtual bool link(const std::string& id, variable_node_i& variable_node) = 0;

        // neighboring variable nodes
        virtual const std::vector<variable_node_i*>& neighbors() const = 0;

protected:
        // check conjugacy of connecting nodes
        virtual bool is_conjugate(size_t i, variable_node_i& variable_node) const = 0;

        // prepare a message to the i'th connected variable
        // node (p messages)
        virtual const p_message_t& message(size_t i) = 0;
};

class variable_node_i : public virtual factor_graph_node_i {
public:
        virtual ~variable_node_i() { }

        virtual variable_node_i* clone() const = 0;

        // get current message
        virtual const q_message_t& operator()() const = 0;

        // send messages to all connected factor nodes
        virtual void send_messages() = 0;

        // link a factor node to this variable node
        virtual mailbox_slot_t<p_message_t>& link(mailbox_slot_t<q_message_t>& slot) = 0;

        // get the type of the distribution this node represents
        virtual const std::type_info& type() const = 0;

protected:
        // prepare the q message
        virtual const q_message_t& message() = 0;
};

// basic implementations of factor and variable nodes
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
                observable_t(factor_node),
        // do not copy the mailer and mailbox, since they should be
        // populated manually to create a new network
                _inbox      (D),
                outbox      (D),
                _neighbors  (D, NULL) {
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
                std::cout << std::endl;
        }
        virtual bool link(size_t i, variable_node_i& variable_node) {
                assert(i < D);
                // allow only conjugate nodes to connect
                if (// variable_node either has to be a conjugate distribution
                    !is_conjugate(i, variable_node) &&
                    // or a dirac distribution
                    variable_node.type() != typeid(dirac_distribution_t)) {
                        return false;
                }
                // pointer to the method that notifies the factor
                // graph about an update
                void (observable_t::*tmp) () const = &factor_node_t::notify;
                // observe the mailbox slot
                _inbox[i].observe(boost::bind(tmp, this));
                // exchange mailbox slots
                outbox[i] = variable_node.link(_inbox[i]);
                // initialize outbox
                outbox[i]->replace(initial_message(i));
                // save neighbor
                _neighbors[i] = &variable_node;
                // return that the nodes were successfully linked
                return true;
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
        void lock_inbox() {
                for (size_t i = 0; i < D; i++) {
                        _inbox[i].lock();
                }
        }
        void unlock_inbox() {
                for (size_t i = 0; i < D; i++) {
                        _inbox[i].unlock();
                }
        }
};

template <typename T>
class variable_node_t : public virtual variable_node_i, public observable_t {
public:
        variable_node_t() {
        }
        variable_node_t(const variable_node_t& variable_node) :
                observable_t   (variable_node),
        // do not copy the inbox and outbox, since they should be
        // populated manually to create a new network
                _inbox         (),
                outbox         (),
                current_message(variable_node.current_message),
                messages       () {
        }

        friend void swap(variable_node_t& left, variable_node_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left._inbox, right._inbox);
                swap(left.outbox, right.outbox);
                swap(left.current_message, right.current_message);
        }

        variable_node_t& operator=(const variable_node_t& node) {
                using std::swap;
                variable_node_t tmp(node);
                swap(*this, tmp);
                return *this;
        }

        virtual const T& operator()() const {
                return current_message();
        }
        virtual void send_messages() {
                std::cout << "variable node " << this << " is sending messages" << std::endl;
                // compute new q-message
                current_message = message();
                // check if this message was sent before
                if (!current_message) {
                        std::cout << "-> new message is the same, stopping." << std::endl;
                        return;
                }
                // lock all connected nodes
                for (size_t i = 0; i < outbox.size(); i++) {
                        std::cout << "-> sending message to neighbor " << i << std::endl;
                        outbox[i]->lock();
                        messages[i] = current_message();
                        // the link is already present
                        outbox[i]->replace(messages[i]);
                        outbox[i]->notify();
                        outbox[i]->unlock();
                }
                std::cout << std::endl;
        }
        virtual mailbox_slot_t<p_message_t>& link(mailbox_slot_t<q_message_t>& slot) {
                size_t i = _inbox.size();
                void (observable_t::*tmp) () const = &variable_node_t::notify;
                // allocate a new message
                messages.push_back(new T());
                // save slot to the outbox
                outbox.push_back(slot);
                // put the current message into the box
                slot.replace(messages[i]);
                // and prepare a new inbox for this node
                _inbox++;
                _inbox[i].observe(boost::bind(tmp, this));
                return _inbox[i];
        }
        virtual const std::type_info& type() const {
                return typeid(T);
        }

protected:
        // prepare the q message
        virtual const T& message() = 0;
        // mailboxes
        _inbox_t<p_message_t> _inbox;
        outbox_t<q_message_t> outbox;
        // messages
        hotnews_t<T> current_message;
        // keep a message for each node
        boost::ptr_vector<T> messages;
};

// specializations of the variable node
////////////////////////////////////////////////////////////////////////////////

template <typename T>
class exponential_vnode_t : public variable_node_t<T> {
public:
        typedef variable_node_t<T> base_t;

        exponential_vnode_t() :
                base_t() {
        }
        exponential_vnode_t(const exponential_vnode_t& exponential_vnode) :
                base_t(exponential_vnode) {
        }
        virtual exponential_vnode_t* clone() const {
                return new exponential_vnode_t(*this);
        }
protected:
        virtual const T& message() {
                // reset new_message
                new_message = T();
                // loop over all slots of the mailbox
                for (size_t i = 0; i < this->_inbox.size(); i++) {
                        // lock this slot
                        this->_inbox[i].lock();
                        // and get the message
                        std::cout << "-> getting message from slot " << i << std::endl;
                        new_message *= this->_inbox[i]();
                        // release lock
                        this->_inbox[i].unlock();
                }
                // normalize message
                new_message.renormalize();
                std::cout << "-> message mean: " << new_message.template moment<1>() << std::endl;
                // has this message be sent before?
                return new_message;
        }
        T new_message;
};

class data_vnode_t : public variable_node_t<dirac_distribution_t> {
public:
        typedef variable_node_t<dirac_distribution_t> base_t;

        data_vnode_t(double x) :
                base_t(),
                data(x) {
        }
        data_vnode_t(const data_vnode_t& data_vnode) :
                base_t(data_vnode),
                data  (data_vnode.data) {
        }
        virtual data_vnode_t* clone() const {
                return new data_vnode_t(*this);
        }
protected:
        virtual const dirac_distribution_t& message() {
                new_message = dirac_distribution_t(data);
                return new_message;
        }
        double data;
        dirac_distribution_t new_message;
};

#endif /* __TFBAYES_FG_NODE_TYPES_HH__ */
