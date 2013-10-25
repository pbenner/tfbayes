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

#ifndef INDEX_HH
#define INDEX_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cstddef>
#include <iostream>

#include <boost/utility.hpp>
#include <boost/functional/hash.hpp> 

#include <tfbayes/dpm/clonable.hh>

// index_t and range_t
////////////////////////////////////////////////////////////////////////////////

class index_i : public clonable, boost::noncopyable {
public:
        virtual index_i* clone() const = 0;
        virtual const size_t& operator[](size_t i) const = 0;
        virtual size_t& operator[](size_t i) = 0;
        virtual bool operator==(const index_i& index) const = 0;
        virtual bool operator!=(const index_i& index) const = 0;
        virtual index_i& operator=(const index_i& index) = 0;
        virtual void operator++(int i) = 0;
        virtual bool operator<(const index_i& index) const = 0;
        virtual size_t hash() const = 0;
};

inline size_t hash_value(const index_i& index) {
        return index.hash();
}
inline index_i* new_clone(const index_i& index)
{
    return index.clone();
}

class index_t : public index_i {
public:
        inline index_t() {
        }
        inline explicit index_t(size_t x0) GCC_ATTRIBUTE_HOT
                : _x0(x0) {
        }
        inline index_t(const index_t& index) GCC_ATTRIBUTE_HOT
                : _x0(index._x0) {
        }

        virtual index_t* clone() const GCC_ATTRIBUTE_HOT {
                return new index_t(*this);
        }

        friend std::ostream& operator<< (std::ostream& o, const index_t& index);

        virtual inline const size_t& operator[](size_t i) const GCC_ATTRIBUTE_HOT {
                return _x0;
        }
        virtual inline size_t& operator[](size_t i) GCC_ATTRIBUTE_HOT {
                return _x0;
        }
        virtual inline bool operator==(const index_i& index) const GCC_ATTRIBUTE_HOT {
                return _x0 == index[0];
        }
        virtual inline bool operator!=(const index_i& index) const GCC_ATTRIBUTE_HOT {
                return _x0 != index[0];
        }
        virtual inline index_t& operator=(const index_i& index) GCC_ATTRIBUTE_HOT {
                _x0 = index[0];
                return *this;
        }
        virtual inline void operator++(int i) {
                _x0++;
        }
        virtual inline bool operator<(const index_i& index) const {
                return _x0 < index[0];
        }
        virtual size_t hash() const {
                boost::hash<size_t> hasher;

                return hasher(_x0);
        }

protected:
        size_t _x0;
};

class seq_index_t : public index_i {
public:
        inline seq_index_t() {
        }
        inline seq_index_t(size_t x0, size_t x1) GCC_ATTRIBUTE_HOT {
                _x[0] = x0;
                _x[1] = x1;
        }
        inline seq_index_t(const seq_index_t& seq_index) GCC_ATTRIBUTE_HOT {
                _x[0] = seq_index[0];
                _x[1] = seq_index[1];
        }

        virtual seq_index_t* clone() const GCC_ATTRIBUTE_HOT {
                return new seq_index_t(*this);
        }

        friend std::ostream& operator<< (std::ostream& o, const seq_index_t& seq_index);

        virtual inline const size_t& operator[](size_t i) const GCC_ATTRIBUTE_HOT {
                return _x[i];
        }
        virtual inline size_t& operator[](size_t i) GCC_ATTRIBUTE_HOT {
                return _x[i];
        }
        virtual inline bool operator==(const index_i& seq_index) const GCC_ATTRIBUTE_HOT {
                return _x[0] == seq_index[0] && _x[1] == seq_index[1];
        }
        virtual inline bool operator!=(const index_i& seq_index) const GCC_ATTRIBUTE_HOT {
                return _x[0] != seq_index[0] || _x[1] != seq_index[1];
        }
        virtual inline seq_index_t& operator=(const index_i& seq_index) GCC_ATTRIBUTE_HOT {
                _x[0] = seq_index[0];
                _x[1] = seq_index[1];
                return *this;
        }
        virtual inline void operator++(int i) {
                _x[1]++;
        }
        virtual inline bool operator<(const index_i& seq_index) const {
                return _x[0] < seq_index[0] || (_x[0] == seq_index[0] && _x[1] < seq_index[1]);
        }
        virtual size_t hash() const {
                size_t seed = 0; 
                boost::hash<size_t> hasher;
                boost::hash_combine(seed, hasher(_x[0]));
                boost::hash_combine(seed, hasher(_x[1]));
                return seed;
        }

protected:
        size_t _x[2];
};

class range_t {
public:
        range_t(const index_i& index, size_t length) __attribute__((hot))
                : _index(index.clone()), _length(length) {
        }
        range_t(const range_t& range) __attribute__((hot))
                : _index (range.index().clone()),
                  _length(range.length()) {
        }
        ~range_t() GCC_ATTRIBUTE_HOT {
                delete(_index);
        }
        range_t& operator=(const range_t& range) GCC_ATTRIBUTE_HOT {
                range_t tmp(range);
                swap(*this, tmp);
                return *this;
        }
        friend void swap(range_t& first, range_t& second) {
                using std::swap;
                swap(first._index,  second._index );
                swap(first._length, second._length);
        }
        bool operator==(const range_t& range) const GCC_ATTRIBUTE_HOT {
                return range.index() == *_index && range.length() == _length;
        }
        const index_i& index() const GCC_ATTRIBUTE_HOT {
                return *_index;
        }
        const size_t& length() const GCC_ATTRIBUTE_HOT {
                return _length;
        }

protected:
        index_i* _index;
        size_t _length;
};

static inline
size_t hash_value(const range_t& range)
{
        return range.index().hash();
}

#endif /* INDEX_HH */
