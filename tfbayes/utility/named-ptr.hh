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

#ifndef __TFBAYES_UTILITY_NAMED_PTR_HH__
#define __TFBAYES_UTILITY_NAMED_PTR_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

template <typename T>
class named_ptr_t : public std::pair<std::string, T*> {
public:
        typedef std::pair<std::string, T*> base_t;

        named_ptr_t() :
                base_t("", NULL)
                { }
        named_ptr_t(std::string tag, T* link) :
                base_t(tag, link)
                { }

        named_ptr_t& operator=(const named_ptr_t& neighbor) {
                base_t::first  = neighbor.first;
                base_t::second = neighbor.second;
                return *this;
        }
        named_ptr_t& operator=(const std::string& lhs) {
                base_t::first = lhs;
                return *this;
        }
        named_ptr_t& operator=(T* lhs) {
                base_t::second = lhs;
                return *this;
        }
        T* operator->() const {
                return base_t::second;
        }
        operator std::string() const {
                return base_t::first;
        }
        operator T*() const {
                return base_t::second;
        }
};

template <typename T>
class named_ptr_vector_t : public std::vector<named_ptr_t<T> > {
public:
        typedef std::vector<named_ptr_t<T> > base_t;

        named_ptr_vector_t() :
                base_t()
                { }
        named_ptr_vector_t(size_t k) :
                base_t(k)
                { }
        named_ptr_vector_t(size_t k, const char *tags[]) :
                base_t(k) {
                for (size_t i = 0; i < k; i++) {
                        base_t::operator[](i) = std::string(tags[i]);
                }
        }
        ssize_t index(const std::string& tag) const {
                for (typename base_t::const_iterator it = base_t::begin();
                     it != base_t::end(); it++) {
                        if (static_cast<std::string>(*it) == tag) {
                                return it - base_t::begin();
                        }
                }
                return -1;
        }
        using base_t::push_back;
        void push_back(const std::string& tag) {
                base_t::push_back(named_ptr_t<T>(tag, NULL));
        }
};

#endif /* __TFBAYES_UTILITY_NAMED_PTR_HH__ */
