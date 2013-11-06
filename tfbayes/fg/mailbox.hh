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

#include <boost/thread.hpp>
#include <boost/function.hpp>
#include <boost/optional.hpp>

template <typename T>
class mailbox_slot_t : public boost::mutex {
public:
        mailbox_slot_t(const boost::function<void ()> f) :
                _fptr_(f) {
        }

        // access the contents
        const T& receive() const {
                return *_ptr_;
        }
        // replace contents
        void replace(const T& ptr) {
                _ptr_ = &ptr;
        }
        // notify the owner of this mailbox slot that a new message
        // has been received
        void notify() const {
                _fptr_();
        }
protected:
        const T* _ptr_;
        boost::function<void ()> _fptr_;
};
template <typename T>
class _inbox_t : public std::vector<boost::optional<mailbox_slot_t<T>&> > {
public:
        typedef std::vector<boost::optional<mailbox_slot_t<T>&> > base_t;

        _inbox_t() :
                base_t() {
        }
        _inbox_t(size_t n) :
                base_t(n, boost::optional<mailbox_slot_t<T>&>()) {
        }
        ~_inbox_t() {
                for (size_t i = 0; i < this->size(); i++) {
                        if (this->operator[](i)) {
                                delete(&*this->operator[](i));
                        }
                }
        }
};
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
