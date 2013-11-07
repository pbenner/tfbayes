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

#ifndef __TFBAYES_FG_HOTNEWS_HH__
#define __TFBAYES_FG_HOTNEWS_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm> /* std::swap */
#include <utility>   /* std::swap */

#include <iostream>

template <typename T>
class hotnews_t : public T {
public:
        hotnews_t() :
                hot(false) {
        }

        friend void swap(hotnews_t& left, hotnews_t& right) {
                using std::swap;
                swap(static_cast<T&>(left), static_cast<T&>(right));
                swap(left.hot, right.hot);
        }

        // standard assignment operator
        hotnews_t& operator=(const hotnews_t& rhs) {
                using std::swap;
                hotnews_t tmp(rhs);
                swap(*this, rhs);
                return *this;
        }
        hotnews_t& operator=(const T& rhs) {
                // if we receive the same news twice, it starts to
                // stink a little
                hot = (*this != rhs);
                if (hot) {
                        std::cout << "---> message is hot" << std::endl;
                        T::operator=(rhs);
                }
                else {
                        std::cout << "---> message is not hot" << std::endl;
                }
                return *this;
        }
        // is this news still hot or does it stink like old fish?
        operator bool() const {
                return hot;
        }
protected:
        bool hot;
};

#endif /* __TFBAYES_FG_HOTNEWS_HH__ */
