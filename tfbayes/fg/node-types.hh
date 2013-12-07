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
#include <boost/lambda/lambda.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <tfbayes/fg/hotnews.hh>
#include <tfbayes/fg/messages.hh>
#include <tfbayes/utility/default-operator.hh>
#include <tfbayes/utility/debug.hh>
#include <tfbayes/utility/linalg.hh>
#include <tfbayes/utility/named-ptr.hh>
#include <tfbayes/utility/observable.hh>

// links between nodes
////////////////////////////////////////////////////////////////////////////////

// access to nodes' messages is implemented by function pointers
typedef boost::function<const p_message_t&()> p_link_t;
typedef boost::function<const q_message_t&()> q_link_t;

// transform message before it is sent to a neighboring node
typedef boost::function<p_message_t&(p_message_t&)> p_map_t;

template <typename T>
class links_t : public std::vector<T> {
public:
        typedef std::vector<T> base_t;

        links_t() :
                base_t()
                { }
        links_t(size_t n) :
                base_t(n, T())
                { }
};

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
// and be able to compute it's contribution to the free energy
class fg_node_i : public virtual clonable, public virtual observable_i {
public:
        virtual ~fg_node_i() { }

        // a node might have an identifier
        virtual const std::string& name() const = 0;

        // execute learning algorithm
        virtual double operator()() = 0;

        // initialize node randomly
        virtual double init(boost::random::mt19937& generator) = 0;
};

class   factor_node_i;
class variable_node_i;

class factor_node_i : public virtual fg_node_i {
public:
        typedef named_ptr_vector_t<variable_node_i> neighbors_t;

        virtual ~factor_node_i() { }

        virtual factor_node_i* clone() const = 0;

        // link a variable node to this factor node
        virtual bool link(const std::string& id, variable_node_i& variable_node,
                          p_map_t f = boost::lambda::_1) = 0;

        // neighboring variable nodes
        virtual const neighbors_t& neighbors() const = 0;

        // notify neighbors
        virtual void notify_neighbors(const variable_node_i& variable_node) const = 0;
};

class variable_node_i : public virtual fg_node_i {
public:
        typedef std::vector<factor_node_i*> neighbors_t;

        virtual ~variable_node_i() { }

        virtual variable_node_i* clone() const = 0;

        // access methods
        virtual const q_message_t& message() const = 0;
        virtual const exponential_family_i& distribution() const = 0;

        // add some data to the variable node
        virtual void condition(const std::matrix<double>& x) = 0;

        // link a factor node to this variable node
        virtual q_link_t link(factor_node_i& factor_node, p_link_t f) = 0;

        // get the type of the distribution this node represents
        virtual const std::type_info& type() const = 0;

        // notify neighbors
        virtual void notify_neighbors() const = 0;

        // neighboring factor nodes
        virtual const neighbors_t& neighbors() const = 0;
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

class factor_node_t : public factor_node_i, public observable_t {
public:
        factor_node_t(size_t k, const char* tags[], const std::string& name = "") :
                _links       (k),
                _neighbors   (k, tags),
                _name        (name) {
                debug(boost::format("allocating factor node %s:%x\n") 
                      % name % this);
        }
        factor_node_t(const factor_node_t& factor_node) :
                factor_node_i(factor_node),
                // do not copy the links, since they should be
                // populated manually to create a new network
                _links       (factor_node._links.size()),
                _neighbors   (factor_node._neighbors),
                _name        (factor_node._name) {
                debug(boost::format("copying factor node %s:%x to %x\n") 
                      % name() % &factor_node % this);
        }
        ~factor_node_t() {
                debug(boost::format("freeing factor node %s:%x\n") 
                      % name() % this);
        }

        virtual factor_node_t* clone() const = 0;

        friend void swap(factor_node_t& left, factor_node_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left._neighbors, right._neighbors);
                swap(left._name,      right._name);
                swap(left._links,     right._links);
        }

        virtual bool link(const std::string& tag, variable_node_i& variable_node, p_map_t f) {
                ssize_t i = _neighbors.index(tag);
                if (i == -1) {
                        return false;
                }
                return link(i, variable_node, f);
        }
        virtual const std::string& name() const {
                return _name;
        }
        virtual void notify_neighbors(const variable_node_i& variable_node) const {
                // notify all other neighbors about an update at node i
                for (size_t i = 0; i < _neighbors.size(); i++) {
                        if (_neighbors[i] != NULL && _neighbors[i] != &variable_node) {
                                _neighbors[i]->notify();
                        }
                }
                observable_t::notify();
        }
        virtual const neighbors_t& neighbors() const {
                return _neighbors;
        }
        virtual double init(boost::random::mt19937& generator) {
                return operator()();
        }
protected:
        bool link(size_t i, variable_node_i& variable_node, p_map_t f) {
                assert(i < _links.size());
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
                // exchange mailbox slots
                _links[i] = variable_node.link(*this, boost::bind(f, boost::bind(&factor_node_t::message, this, i)));
                // save neighbor
                _neighbors[i] = &variable_node;
                // return that the nodes were successfully linked
                return true;
        }
        // prepare a message to the i'th connected variable
        // node (p messages)
        virtual p_message_t& message(size_t i) = 0;
        // check conjugacy of connecting nodes
        virtual bool is_conjugate(size_t i, variable_node_i& variable_node) const = 0;
        // links to neighboring nodes
        links_t<q_link_t> _links;
        // keep track of neighboring nodes for cloning whole networks
        neighbors_t _neighbors;
        // id of this node
        std::string _name;
};

template <typename T>
class variable_node_t : public variable_node_i, public observable_t {
public:
        variable_node_t(const T& distribution, const std::string& name = "") :
                variable_node_i(),
                observable_t   (),
                _name          (name),
                _distribution  (distribution),
                _message       (q_message_t(distribution.k())) {
                debug(boost::format("allocating variable node %s:%x\n") 
                      % name % this);
                if (distribution) {
                        _message = distribution.moments();
                }
        }
        variable_node_t(const variable_node_t& variable_node) :
                variable_node_i(variable_node),
                observable_t   (variable_node),
                _links         (),
                _neighbors     (),
                _name          (variable_node._name),
                _distribution  (variable_node._distribution),
                _message       (variable_node._message) {
                debug(boost::format("copying variable node %s:%x to %x\n") 
                      % name() % &variable_node % this);
        }
        ~variable_node_t() {
                debug(boost::format("freeing variable node %s:%x\n") 
                      % name() % this);
        }
        // only derived classes should be cloned
        virtual variable_node_t* clone() const {
                return new variable_node_t(*this);
        }

        friend void swap(variable_node_t& left, variable_node_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left._name,         right._name);
                swap(left._links,        right._links);
                swap(left._neighbors,    right._neighbors);
                swap(left._distribution, right._distribution);
                swap(left._message,      right._message);
        }
        virtual_assignment_operator(variable_node_t)
        derived_assignment_operator(variable_node_t, variable_node_i)

        virtual q_link_t link(factor_node_i& factor_node, p_link_t f) {
                _links.push_back(f);
                _neighbors.push_back(&factor_node);
                // return a lambda to the call operator
                return boost::bind(&variable_node_t::message, this);
        }
        virtual const std::type_info& type() const {
                return typeid(T);
        }
        virtual const std::string& name() const {
                return _name;
        }
        virtual const neighbors_t& neighbors() const {
                return _neighbors;
        }
        virtual const q_message_t& message() const {
                return static_cast<const q_message_t&>(_message);
        }
        virtual double operator()() {
                // compute new q-message
                debug(boost::format("variable node %s:%x is preparing a new message\n")
                      % name() % this);
                // update distribution
                _distribution.clear();
                for (size_t i = 0; i < _links.size(); i++) {
                        // get the message
                        _distribution *= _links[i]();
                }
                // normalize message
                _distribution.renormalize();
                // the new message is the moments of the
                // sufficient statistics
                _message = _distribution.moments();
                debug(boost::format("variable node %s:%x has new message: ")
                      % name() % this << std::boolalpha
                      << (bool)_message << std::endl);
                // check if this message was sent before
                if (_message) notify_neighbors();
                debug(boost::format("variable node %s:%x computed free energy: %d\n")
                      % name() % this % _distribution.entropy());
                return _distribution.entropy();
        }
        virtual void condition(const std::matrix<double>& x);
        virtual const T& distribution() const {
                return _distribution;
        }
        void notify_neighbors() const {
                for (size_t i = 0; i < neighbors().size(); i++) {
                        neighbors()[i]->notify_neighbors(*this);
                }
        }
        virtual double init(const T& distribution) {
                _distribution = distribution;
                _message      = distribution.moments();
                return distribution.entropy();
        }
        virtual double init(boost::random::mt19937& generator) {
                return 0.0;
        }
protected:
        // links to neighboring nodes
        links_t<p_link_t> _links;
        // keep track of neighboring nodes for cloning whole networks
        neighbors_t _neighbors;
        // id of this node
        std::string _name;
        // save current distribution to compute the entropy
        T _distribution;
        // current message
        hotnews_t<q_message_t> _message;
};

// variable node specialization for holding data
////////////////////////////////////////////////////////////////////////////////

template <typename T>
class data_vnode_t : public variable_node_t<T> {
public:
        typedef variable_node_t<T> base_t;

        data_vnode_t(const T& distribution, const std::string& name) :
                base_t(distribution, name) {
        }
        data_vnode_t(const variable_node_t<T>& variable_node) :
                base_t(variable_node) {
        }

        virtual data_vnode_t* clone() const {
                return new data_vnode_t(*this);
        }

        virtual double operator()() {
                return 0.0;
        }
        virtual void condition(const std::matrix<double>& x) {
                debug(boost::format("data_vnode %s:%x is receiving new data")
                      % base_t::name() % this << std::endl);
                // reset message
                base_t::_message.clear();
                // loop over data points
                for (size_t i = 0; i < x.size(); i++) {
                        // compute statistics of a single data point
                        statistics_t statistics = base_t::_distribution.statistics(x[i]);
                        // make sure we're talking the same language
                        assert(statistics.size() == base_t::_message.size());
                        // add to current message
                        base_t::_message += statistics;
                }
                assert(x.size() == base_t::_message.n);
        }
        virtual void notify() const {
        }
};

// convert a variable node to a data node
////////////////////////////////////////////////////////////////////////////////

template <typename T>
void
variable_node_t<T>::condition(const std::matrix<double>& x) {
        assert(sizeof(variable_node_t<T>) == sizeof(data_vnode_t<T>));
        // smalltalk "become"
        variable_node_t<T> tmp(*this);
        // save neighbors and links
        swap(tmp, *this);
        this->~variable_node_t<T>();
        new (this) data_vnode_t<T>(tmp);
        swap(tmp, *this);
        this->condition(x);
}

#endif /* __TFBAYES_FG_NODE_TYPES_HH__ */
