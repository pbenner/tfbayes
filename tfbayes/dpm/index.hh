/* Copyright (C) 2011-2013 Philipp Benner
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

#ifndef __TFBAYES_DPM_INDEX_HH__
#define __TFBAYES_DPM_INDEX_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cstddef>
#include <iostream>

#include <boost/utility.hpp>
#include <boost/functional/hash.hpp> 

// index_t and range_t
////////////////////////////////////////////////////////////////////////////////


class index_t {
public:
        inline index_t() {
        }
        inline index_t(ssize_t x0) GCC_ATTRIBUTE_HOT {
                m_x[0] = x0;
                m_x[1] = -1;
        }
        inline index_t(ssize_t x0, ssize_t x1) GCC_ATTRIBUTE_HOT {
                m_x[0] = x0;
                m_x[1] = x1;
        }
        inline index_t(const index_t& index) GCC_ATTRIBUTE_HOT {
                m_x[0] = index[0];
                m_x[1] = index[1];
        }

        friend std::ostream& operator<< (std::ostream& o, const index_t& index);

        inline const ssize_t& operator[](size_t i) const GCC_ATTRIBUTE_HOT {
                return m_x[i];
        }
        inline ssize_t& operator[](size_t i) GCC_ATTRIBUTE_HOT {
                return m_x[i];
        }
        inline void operator++(int i) {
                if (m_x[1] == -1) m_x[0]++;
                else              m_x[1]++;
        }
        inline bool operator==(const index_t& index) const GCC_ATTRIBUTE_HOT {
                return m_x[0] == index[0] && m_x[1] == index[1];
        }
        inline bool operator!=(const index_t& index) const GCC_ATTRIBUTE_HOT {
                return m_x[0] != index[0] || m_x[1] != index[1];
        }
        inline bool operator<(const index_t& index) const {
                return m_x[0] < index[0] || (m_x[0] == index[0] && m_x[1] < index[1]);
        }
        inline size_t dim() const {
                if (m_x[1] == -1) return 1;
                else              return 2;
        }
        size_t hash() const {
                size_t seed = 0; 
                boost::hash<size_t> hasher;
                boost::hash_combine(seed, hasher(m_x[0]));
                boost::hash_combine(seed, hasher(m_x[1]));
                return seed;
        }

protected:
        ssize_t m_x[2];
};

class range_t {
public:
        range_t()
                : m_index  ()
                , m_length ()
                , m_reverse()
                { }
        range_t(const index_t& index, size_t length, bool reverse) GCC_ATTRIBUTE_HOT
                : m_index  (index)
                , m_length (length)
                , m_reverse(reverse) {
        }
        range_t(const index_t& index, size_t length) GCC_ATTRIBUTE_HOT
                : m_index  (index)
                , m_length (length)
                , m_reverse(false) {
        }

        friend std::ostream& operator<< (std::ostream& o, const range_t& range);

        bool operator==(const range_t& range) const GCC_ATTRIBUTE_HOT {
                return range.index  () == m_index  &&
                       range.length () == m_length &&
                       range.reverse() == range.reverse();
        }
        bool operator!=(const range_t& range) const GCC_ATTRIBUTE_HOT {
                return !operator==(range);
        }
        const index_t& index() const GCC_ATTRIBUTE_HOT {
                return m_index;
        }
              index_t& index() GCC_ATTRIBUTE_HOT {
                return m_index;
        }
        const size_t& length() const GCC_ATTRIBUTE_HOT {
                return m_length;
        }
              size_t& length() GCC_ATTRIBUTE_HOT {
                return m_length;
        }
        const bool& reverse() const GCC_ATTRIBUTE_HOT {
                return m_reverse;
        }
              bool& reverse() GCC_ATTRIBUTE_HOT {
                return m_reverse;
        }
protected:
        index_t m_index;
        size_t  m_length;
        bool    m_reverse;
};

inline size_t hash_value(const index_t& index) {
        return index.hash();
}

inline size_t hash_value(const range_t& range)
{
        return range.index().hash();
}

#endif /* __TFBAYES_DPM_INDEX_HH__ */
