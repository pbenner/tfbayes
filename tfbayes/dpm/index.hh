/* Copyright (C) 2011 Philipp Benner
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

#include <stddef.h>

#include <clonable.hh>
#include <iostream>

#include <boost/functional/hash.hpp> 

// index_t and range_t
////////////////////////////////////////////////////////////////////////////////

class index_i : public clonable {
public:
        virtual index_i* clone() const = 0;
        virtual const size_t& operator[](size_t i) const = 0;
        virtual size_t& operator[](size_t i) = 0;
        virtual bool operator==(const index_i& index) const = 0;
        virtual index_i& operator=(const index_i& index) = 0;
        virtual void operator++(int i) = 0;
        virtual bool operator<(const index_i& index) const = 0;
        virtual size_t hash() const = 0;
};

class index_t : public index_i {
public:
        __inline__ index_t() {
        }
        __inline__ explicit index_t(size_t x0) : _x0(x0) {
        }
        __inline__ index_t(const index_t& index) : _x0(index._x0) {
        }

        virtual index_i* clone() const { return new index_t(*this); }

        friend std::ostream& operator<< (std::ostream& o, const index_t& index);

        virtual __inline__ const size_t& operator[](size_t i) const {
                return _x0;
        }
        virtual __inline__ size_t& operator[](size_t i) {
                return _x0;
        }
        virtual __inline__ bool operator==(const index_i& index) const {
                return _x0 == index[0];
        }
        virtual __inline__ index_t& operator=(const index_i& index) {
                _x0 = index[0];
                return *this;
        }
        virtual __inline__ void operator++(int i) {
                _x0++;
        }
        virtual __inline__ bool operator<(const index_i& index) const {
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
        __inline__ seq_index_t() {
        }
        __inline__ seq_index_t(size_t x0, size_t x1) {
                _x[0] = x0;
                _x[1] = x1;
        }
        __inline__ seq_index_t(const seq_index_t& seq_index) {
                _x[0] = seq_index[0];
                _x[1] = seq_index[1];
        }

        virtual index_i* clone() const { return new seq_index_t(*this); }

        friend std::ostream& operator<< (std::ostream& o, const seq_index_t& seq_index);

        virtual __inline__ const size_t& operator[](size_t i) const {
                return _x[i];
        }
        virtual __inline__ size_t& operator[](size_t i) {
                return _x[i];
        }
        virtual __inline__ bool operator==(const index_i& seq_index) const {
                return _x[0] == seq_index[0] && _x[1] == seq_index[1];
        }
        virtual __inline__ seq_index_t& operator=(const index_i& seq_index) {
                _x[0] = seq_index[0];
                _x[1] = seq_index[1];
                return *this;
        }
        virtual __inline__ void operator++(int i) {
                _x[1]++;
        }
        virtual __inline__ bool operator<(const index_i& seq_index) const {
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
        range_t(const index_i& index, size_t length) 
                : _index(index.clone()), _length(length) {
        }
        range_t(const range_t& range) {
                operator=(range);
        }
        ~range_t() {
                delete(_index);
        }
        range_t& operator=(const range_t& range) {
                _index  = range.index().clone();
                _length = range.length();
                return *this;
        }
        bool operator==(const range_t& range) const {
                return range.index() == *_index && range.length() == _length;
        }
        const index_i& index() const {
                return *_index;
        }
        const size_t& length() const {
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
