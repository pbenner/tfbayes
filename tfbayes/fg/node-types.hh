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

#include <tfbayes/fg/hotnews.hh>
#include <tfbayes/fg/mailbox.hh>
#include <tfbayes/fg/messages.hh>
#include <tfbayes/utility/default-operator.hh>
#include <tfbayes/utility/debug.hh>
#include <tfbayes/utility/observable.hh>

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
class fg_node_i : public virtual clonable, public virtual observable_i {
public:
        virtual ~fg_node_i() { }

        // send messages to all connected variable nodes
        virtual void send_messages() = 0;

        // compute free energy
        virtual double free_energy() const = 0;

        // a node might have an identifier
        virtual const std::string& name() const = 0;
};

class   factor_node_i;
class variable_node_i;

class factor_node_i : public fg_node_i {
public:
        virtual ~factor_node_i() { }

        virtual factor_node_i* clone() const = 0;

        virtual factor_node_i& operator=(const factor_node_i& factor_node) = 0;

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

class variable_node_i : public virtual fg_node_i {
public:
        virtual ~variable_node_i() { }

        virtual variable_node_i* clone() const = 0;

        virtual variable_node_i& operator=(const variable_node_i& variable_node) = 0;

        // get current message
        virtual const exponential_family_i& operator()() const = 0;

        // link a factor node to this variable node
        virtual mailbox_slot_t<p_message_t>& link(mailbox_slot_t<q_message_t>& slot) = 0;

        // get the type of the distribution this node represents
        virtual const std::type_info& type() const = 0;

        // condition this node on some data
        virtual void condition(const std::vector<double>& x) = 0;

protected:
        // prepare the q message
        virtual const q_message_t& message() = 0;
};

// tell boost how to clone nodes
inline factor_node_i* new_clone(const factor_node_i& a)
{
    return a.clone();
}
inline variable_node_i* new_clone(const variable_node_i& a)
{
    return a.clone();
}

// basic implementations of factor and variable nodes
////////////////////////////////////////////////////////////////////////////////

template <size_t D>
class factor_node_t : public factor_node_i, public observable_t {
public:
        factor_node_t(const std::string& name = "") :
                _inbox       (D),
                outbox       (D),
                _neighbors   (D, NULL),
                _name        (name) {
                debug("allocating factor node at " << this << std::endl);
        }
        factor_node_t(const factor_node_t& factor_node) :
                factor_node_i(factor_node),
                observable_t (factor_node),
                // do not copy the mailer and mailbox, since they should be
                // populated manually to create a new network
                _inbox       (factor_node._inbox),
                outbox       (D),
                _neighbors   (D, NULL),
                _name        (factor_node._name) {
                debug("copying factor node from " << &factor_node << " to " << this << std::endl);
        }

        virtual factor_node_t* clone() const = 0;

        friend void swap(factor_node_t& left, factor_node_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left._neighbors, right._neighbors);
                swap(left._name,      right._name);
                swap(left._inbox,     right._inbox);
        }
        factor_node_t& operator=(const factor_node_t& factor_node) = delete;

        virtual void send_messages() {
                for (size_t i = 0; i < D; i++) {
                        if (outbox[i]) {
                                debug(boost::format("factor node %s:%x is sending a message "
                                                    "to variable node %s:%x\n")
                                      % name() % this % _neighbors[i]->name() % _neighbors[i]);
                                (*outbox[i])() = message(i);
                                  outbox[i]->notify();
                                debug("................................................................................"
                                      << std::endl);
                        }
                }
        }
        virtual bool link(size_t i, variable_node_i& variable_node) {
                assert(i < D);
                // allow only conjugate nodes to connect
                debug(boost::format("attempting to link factor node %s:%x "
                                    " with variable node %s:%x")
                      % name() % this % variable_node.name() % &variable_node
                      << std::endl);
                if (// variable_node has to be a conjugate distribution
                    !is_conjugate(i, variable_node)) {
                        debug("-> failed!" << std::endl);
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
                (*outbox[i])() = initial_message(i);
                // save neighbor
                _neighbors[i] = &variable_node;
                // return that the nodes were successfully linked
                return true;
        }
        virtual const std::vector<variable_node_i*>& neighbors() const {
                return _neighbors;
        }
        virtual const std::string& name() const {
                return _name;
        }
protected:
        // initialize outbox i
        virtual const p_message_t& initial_message(size_t i) const = 0;
        // mailboxes
        _inbox_t<q_message_t> _inbox;
        outbox_t<p_message_t> outbox;
        // keep track of neighboring nodes for cloning whole networks
        std::vector<variable_node_i*> _neighbors;
        // id of this node
        std::string _name;
};

template <typename T>
class variable_node_t : public variable_node_i, public observable_t {
public:
        variable_node_t(const std::string& name = "") :
                _name          (name) {
                debug("allocating variable node at " << this << std::endl);
        }
        variable_node_t(const variable_node_t& variable_node) :
                variable_node_i(variable_node),
                observable_t   (variable_node),
        // do not copy the inbox and outbox, since they should be
        // populated manually to create a new network
                _inbox         (),
                outbox         (),
                current_message(variable_node.current_message),
                _name          (variable_node._name) {
                debug("copying variable node from " << &variable_node << " to " << this << std::endl);
        }

        friend void swap(variable_node_t& left, variable_node_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left.current_message, right.current_message);
                swap(left._name,           right._name);
        }
        variable_node_t& operator=(const variable_node_t& variable_node) = delete;

        virtual void send_messages() {
                // compute new q-message
                debug(boost::format("variable node %s:%x is preparing a new message\n")
                      % name() % this);
                current_message = message();
                // check if this message was sent before
                if (!current_message) {
                        debug(boost::format("variable node %s:%x has no new message\n")
                              % name() % this);
                        debug("--------------------------------------------------------------------------------"
                              << std::endl);
                        return;
                }
                debug(boost::format("variable node %s:%x is sending messages\n")
                      % name() % this);
                for (size_t i = 0; i < outbox.size(); i++) {
                        (*outbox[i])() = current_message;
                          outbox[i]->notify();
                }
        }
        virtual mailbox_slot_t<p_message_t>& link(mailbox_slot_t<q_message_t>& slot) {
                size_t i = _inbox.size();
                void (observable_t::*tmp) () const = &variable_node_t::notify;
                // save slot to the outbox
                outbox.push_back(slot);
                // and prepare a new inbox for this node
                _inbox++;
                _inbox[i].replace(new T());
                _inbox[i].observe(boost::bind(tmp, this));
                return _inbox[i];
        }
        virtual const std::type_info& type() const {
                return typeid(T);
        }
        virtual const std::string& name() const {
                return _name;
        }
protected:
        // prepare the q message
        virtual const typename T::moments_t& message() = 0;
        // mailboxes
        _inbox_t<p_message_t> _inbox;
        outbox_t<q_message_t> outbox;
        // messages
        hotnews_t<typename T::moments_t> current_message;
        // id of this node
        std::string _name;
};

// specializations of the variable node
////////////////////////////////////////////////////////////////////////////////

template <typename T>
class exponential_vnode_t : public variable_node_t<T> {
public:
        typedef variable_node_t<T> base_t;

        exponential_vnode_t(const std::string& name = "") :
                base_t(name) {
        }
        exponential_vnode_t(const exponential_vnode_t& exponential_vnode) :
                base_t      (exponential_vnode),
                new_message (exponential_vnode.new_message) {
        }
        virtual exponential_vnode_t* clone() const {
                return new exponential_vnode_t(*this);
        }
        friend void swap(exponential_vnode_t& left, exponential_vnode_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
                swap(left.new_message, right.new_message);
        }
        virtual_assignment_operator(exponential_vnode_t)
        derived_assignment_operator(exponential_vnode_t, variable_node_i)

        virtual void condition(const std::vector<double>& x) {
        }
        virtual const T& operator()() const {
                return distribution;
        }
        virtual double free_energy() const {
                debug(boost::format("variable node %s:%x computed entropy: %d\n")
                      % base_t::name() % this % distribution.entropy());
                debug("--------------------------------------------------------------------------------"
                      << std::endl);
                return distribution.entropy();
        }
protected:
        virtual const typename T::moments_t& message() {
                // get a new exponential family
                distribution = T();
                // loop over all slots of the mailbox
                for (size_t i = 0; i < this->_inbox.size(); i++) {
                        // get the message
                        distribution *= this->_inbox[i]();
                }
                // normalize message
                distribution.renormalize();
                debug(std::endl);
                // the new message is the moments of the exponential
                // family
                new_message = distribution.moments();
                // has this message be sent before?
                return new_message;
        }
        T distribution;
        typename T::moments_t new_message;
};

template <typename T>
class data_vnode_t : public variable_node_t<T> {
public:
        typedef variable_node_t<T> base_t;

        data_vnode_t(const std::string& name) :
                base_t(name),
                new_message() {
        }
        data_vnode_t(const data_vnode_t& data_vnode) :
                base_t      (data_vnode),
                new_message (data_vnode.new_message) {
        }
        virtual data_vnode_t* clone() const {
                return new data_vnode_t(*this);
        }

        friend void swap(data_vnode_t& left, data_vnode_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
                swap(left.new_message, right.new_message);
        }
        virtual_assignment_operator(data_vnode_t)
        derived_assignment_operator(data_vnode_t, variable_node_i)

        void condition(const std::vector<double>& x) {
                debug(boost::format("data_vnode %s:%x is receiving new data")
                      % base_t::name() % this << std::endl);
                new_message = dirac_distribution_t(x);
        }
        virtual const T& operator()() const {
                return distribution;
        }
        virtual double free_energy() const {
                debug("--------------------------------------------------------------------------------"
                      << std::endl);
                return 0.0;
        }
protected:
        virtual const typename T::moments_t& message() {
                return new_message;
        }
        // this is only a dummy distribution that is not needed
        T distribution;
        // the actual message
        typename T::moments_t new_message;
};

#endif /* __TFBAYES_FG_NODE_TYPES_HH__ */
