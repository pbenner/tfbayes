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

#ifndef __TFBAYES_UTILITY_OBSERVABLE_HH__
#define __TFBAYES_UTILITY_OBSERVABLE_HH__

#include <algorithm> /* std::swap */
#include <utility>   /* std::swap */

class observable_i {
public:
        virtual ~observable_i() { }

        virtual void notify() const = 0;
        virtual void observe(const boost::function<void ()>& f) = 0;
};

class observable_t : public virtual observable_i {
public:
        virtual ~observable_t() { }

        friend void swap(observable_t& left, observable_t& right) {
                using std::swap;
                swap(left._notify, right._notify);
        }

        virtual void notify() const {
                if (_notify) _notify();
        }
        virtual void observe(const boost::function<void ()>& f) {
                _notify = f;
        }
private:
        boost::function<void ()> _notify;
};

#endif /* __TFBAYES_UTILITY_OBSERVABLE_HH__ */
