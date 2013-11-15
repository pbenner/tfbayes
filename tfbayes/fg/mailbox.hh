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

#ifndef __TFBAYES_FG_MAILBOX_HH__
#define __TFBAYES_FG_MAILBOX_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/observable.hh>

template <typename T>
class mailbox_slot_t : public virtual clonable, public observable_t {
public:
        mailbox_slot_t() :
                observable_t(),
                _ptr_       (NULL) {
        }
        mailbox_slot_t(const mailbox_slot_t& mailbox_slot) :
                observable_t(mailbox_slot),
                _ptr_       (NULL) {
                if (mailbox_slot._ptr_ != NULL) {
                        _ptr_ = mailbox_slot._ptr_->clone();
                }
        }
        virtual ~mailbox_slot_t() {
                if (_ptr_ != NULL) {
                        delete(_ptr_);
                }
        }

        mailbox_slot_t* clone() const {
                return new mailbox_slot_t(*this);
        }
        friend void swap(mailbox_slot_t& left, mailbox_slot_t& right) {
                using std::swap;
                swap(static_cast<observable_t&>(left),
                     static_cast<observable_t&>(right));
                swap(left._ptr_, right._ptr_);
        }
        virtual_assignment_operator(mailbox_slot_t)

        // access the contents
        const T& operator()() const {
                return *_ptr_;
        }
        T& operator()() {
                return *_ptr_;
        }
        // replace contents
        void replace(T& ptr) {
                if (_ptr_ != NULL) {
                        delete(_ptr_);
                }
                _ptr_ = &ptr;
        }
        void replace(T* ptr) {
                if (_ptr_ != NULL) {
                        delete(_ptr_);
                }
                _ptr_ = ptr;
        }
protected:
        T* _ptr_;
};

template<typename T>
mailbox_slot_t<T>* new_clone(const mailbox_slot_t<T>& a)
{
    return a.clone();
}

// For variable nodes, it is necessary that the inbox is able to have
// as many slots as factor nodes are connected to the variable nodes,
// however, this might vary and is not determined by the type of the
// variable node. Factor nodes on the other hand, have a fixed number
// of slots that depends on the type of the node. If a variable node
// is connected to a slot, it receives a link to this slot which is
// added to the outbox of the variable node. If no variable node is
// connected to the slot, then a dirac distribution is placed in the
// slot that represents the respective fixed parameter.
template <typename T>
class _inbox_t : public boost::ptr_vector<mailbox_slot_t<T> > {
public:
        typedef boost::ptr_vector<mailbox_slot_t<T> > base_t;

        _inbox_t() :
                base_t() {
        }
        _inbox_t(size_t n) :
                base_t() {
                for (size_t i = 0; i < n; i++)
                        (*this)++;
        }
        // add an additional slot
        _inbox_t& operator++(int) {
                this->push_back(new mailbox_slot_t<T>());
                return *this;
        }
};

// An outbox has a fixed number of slots for factor nodes. If a
// variable node connects to the factor, the respective slot of the
// outbox points to the desired mailbox slot of the connecting node.
// Variable nodes push as many slots to the outbox as needed, the
// number is not determined by the type of the variable node.
template <typename T>
class outbox_t : public std::vector<boost::optional<mailbox_slot_t<T>&> > {
public:
        typedef std::vector<boost::optional<mailbox_slot_t<T>&> > base_t;

        outbox_t() :
                base_t() {
        }
        outbox_t(size_t n) :
                base_t(n, boost::optional<mailbox_slot_t<T>&>()) {
        }
};

#endif /* __TFBAYES_FG_MAILBOX_HH__ */
