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
private:
        // receive a message from a variable node (q message), this
        // method should do nothing but to save the message and notify
        // the factor graph
        // this function must be private so it can only be used by the link
        // function of the class that implements linking nodes
        virtual void recv_message(size_t i, const q_message_t& msg) = 0;
};

class variable_node_i : public virtual factor_graph_node_i {
public:
        virtual ~variable_node_i() { }

        virtual variable_node_i* clone() const = 0;

        // send messages to all connected factor nodes
        virtual void send_messages() = 0;

        // link a factor node to this variable node
        virtual boost::function<void (const p_message_t&)>
        link(boost::function<void (const q_message_t&)> f, boost::mutex& m) = 0;

        // get the lock of this node
        virtual boost::mutex& lock() = 0;

protected:
        // prepare the message to be sent
        virtual const p_message_t& message() = 0;
private:
        // receive a message from a factor node (p message), this
        // method should do nothing but to save the message and notify
        // the factor graph
        // this function must be private so it can only be used by the link
        // function of the class that implements linking nodes
        virtual void recv_message(size_t i, const p_message_t& msg) = 0;
};

// basic implementations
////////////////////////////////////////////////////////////////////////////////

template <size_t D>
class factor_node_t : public virtual factor_node_i, public observable_t {
public:
        factor_node_t() :
                mailer    (),
                mailbox   (),
                _neighbors(D, NULL) {
                // initialize mailbox
                std::fill(mailbox.begin(), mailbox.end(), (const q_message_t*)NULL);
        }
        factor_node_t(const factor_node_t& factor_node) :
        // do not copy the mailer and mailbox, since they should be
        // populated manually to create a new network
                mailer    (),
                mailbox   (),
                _neighbors(D, NULL) {
                // initialize mailbox
                std::fill(mailbox.begin(), mailbox.end(), (const q_message_t*)NULL);
        }

        virtual factor_node_t* clone() const = 0;

        friend void swap(factor_node_t& left, factor_node_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left.mailer,     right.mailer);
                swap(left.mailbox,    right.mailbox);
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
                lock_mailbox();
                for (size_t i = 0; i < D; i++) {
                        if (mailer[i]) {
                                message_locks[i].lock();
                                mailer[i](message(i));
                                message_locks[i].unlock();
                        }
                }
                unlock_mailbox();
        }
        virtual void link(size_t i, variable_node_i& variable_node) {
                assert(i < D);
                // pointer to the method that receives messages
                void (factor_node_t::*tmp) (size_t, const q_message_t&) = &factor_node_t::recv_message;
                // receive and send mailer
                mailer[i] = variable_node.link(
                        boost::bind(tmp, this, i, _1), message_locks[i]);
                // receive the lock of the variable node
                mailbox_locks[i] = variable_node.lock();
                // save neighbor
                _neighbors[i] = &variable_node;
        }
        virtual const std::vector<variable_node_i*>& neighbors() const {
                return _neighbors;
        }

protected:
        boost::array<boost::function<void (const p_message_t&)>, D> mailer;
        boost::array<const q_message_t*, D> mailbox;
        // keep track of neighboring nodes for cloning whole networks
        std::vector<variable_node_i*> _neighbors;
private:
        virtual void recv_message(size_t i, const q_message_t& msg) {
                std::cout << "factor node has received a message from neighbor " << i << std::endl;
                // received a message from neighbor i
                mailbox_locks[i]->lock();
                mailbox[i] = &msg;
                mailbox_locks[i]->unlock();
        }
        // lock every slot in the mailbox, so that we can prepare a
        // new message without receiving new mail while doing so
        void lock_mailbox() const {
                for (size_t i = 0; i < D; i++) {
                        if (mailbox_locks[i]) mailbox_locks[i]->lock();
                }
        }
        void unlock_mailbox() const {
                for (size_t i = 0; i < D; i++) {
                        if (mailbox_locks[i]) mailbox_locks[i]->unlock();
                }
        }
        // locking mechanism for locally saved messages
        mutable boost::array<boost::mutex, D> message_locks;
        // locks for all slots in the mailbox
        mutable boost::array<boost::optional<boost::mutex&>, D> mailbox_locks;
};

template <typename T>
class variable_node_t : public virtual variable_node_i, public observable_t {
public:
        variable_node_t() :
                mailer      (),
                mailbox     (),
                old_message (),
                new_message () {
                // initialize mailbox
                std::fill(mailbox.begin(), mailbox.end(), (const p_message_t*)NULL);
        }
        variable_node_t(const variable_node_t& variable_node) :
        // do not copy the mailer and mailbox, since they should be
        // populated manually to create a new network
                mailer      (),
                mailbox     (),
                old_message (variable_node.old_message),
                new_message (variable_node.new_message) {
                // initialize mailbox
                std::fill(mailbox.begin(), mailbox.end(), (const p_message_t*)NULL);
        }
        virtual variable_node_t* clone() const {
                return new variable_node_t(*this);
        }

        friend void swap(variable_node_t& left, variable_node_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left.mailer,      right.mailer);
                swap(left.mailbox,     right.mailbox);
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
                // check if there are any nodes connected
                if (mailbox.size() == 0) return;
                // prepare the message
                const p_message_t& msg = message();
                // check if this message was sent before
                if (msg == old_message) return;
                // send message
                for (size_t i = 0; i < mailer.size(); i++) {
                        std::cout << "-> sending message to neighbor " << i << std::endl;
                        mailer[i](msg);
                }
                // save message
                old_message = msg;
        }
        virtual boost::function<void (const p_message_t&)> link(boost::function<void (const q_message_t&)> f, boost::mutex& m) {
                void (variable_node_t::*tmp) (size_t, const p_message_t&) = &variable_node_t::recv_message;
                size_t i = mailer.size();
                mailer .push_back(f);
                mailbox.push_back(NULL);
                mailbox_locks.push_back(m);
                return boost::bind(tmp, this, i, _1);
        }
        virtual boost::mutex& lock() {
                return message_lock;
        }

protected:
        virtual const p_message_t& message() {
                // a temporary message
                T msg;
                // loop over all slots of the mailbox
                for (size_t i = 0; i < mailbox.size(); i++) {
                        // is there something in slot i to receive?
                        if (!mailbox_locks[i]) {
                                continue;
                        }
                        // lock this slot
                        mailbox_locks[i]->lock();
                        // and get the message
                        if (mailbox[i]) {
                                msg *= static_cast<const T&>(*mailbox[i]);
                        }
                        // release lock
                        mailbox_locks[i]->unlock();
                }
                // update locally saved message
                message_lock.lock();
                new_message = msg;
                message_lock.unlock();

                return new_message;
        }

        std::vector<boost::function<void (const q_message_t&)> > mailer;
        std::vector<const p_message_t*> mailbox;

        T old_message;
        T new_message;
private:
        virtual void recv_message(size_t i, const p_message_t& msg) {
                std::cout << "variable node has received a message from neighbor " << i << std::endl;
                // received a message from neighbor i
                mailbox_locks[i]->lock();
                mailbox[i] = &msg;
                mailbox_locks[i]->unlock();
        }
        // locking mechanism for the locally saved message
        mutable boost::mutex message_lock;
        // locks for all slots in the mailbox
        mutable std::vector<boost::optional<boost::mutex&> > mailbox_locks;
};

#endif /* __TFBAYES_FG_NODE_TYPES_HH__ */
