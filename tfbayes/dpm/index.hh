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
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stddef.h>

#include <clonable.hh>
#include <iostream>

// index_t and range_t
////////////////////////////////////////////////////////////////////////////////

class index_t : public clonable {
public:
        __inline__ index_t() {
        }
        __inline__ explicit index_t(size_t x0) : _x0(x0) {
        }
        __inline__ index_t(const index_t& index) : _x0(index._x0) {
        }

        index_t* clone() const { return new index_t(*this); }

        friend std::ostream& operator<< (std::ostream& o, const index_t& index);

        virtual __inline__ const size_t& operator[](size_t i) const {
                return _x0;
        }
        virtual __inline__ size_t& operator[](size_t i) {
                return _x0;
        }
        virtual __inline__ bool operator==(const index_t& index) const {
                return _x0 == index[0];
        }
        virtual __inline__ void operator=(const index_t& index) {
                _x0 = index[0];
        }
        virtual __inline__ void operator++(int i) {
                _x0++;
        }

protected:
        size_t _x0;
};

class seq_index_t : public index_t {
public:
        __inline__ seq_index_t() : index_t() {
        }
        __inline__ seq_index_t(size_t x0, size_t x1) : index_t(x0), _x1(x1) {
        }
        __inline__ seq_index_t(const seq_index_t& seq_index) : index_t(seq_index), _x1(seq_index._x1) {
        }

        index_t* clone() const { return new seq_index_t(*this); }

        friend std::ostream& operator<< (std::ostream& o, const seq_index_t& seq_index);

        virtual __inline__ const size_t& operator[](size_t i) const {
                return i == 0 ? _x0 : _x1;
        }
        virtual __inline__ size_t& operator[](size_t i) {
                return i == 0 ? _x0 : _x1;
        }
        virtual __inline__ bool operator==(const seq_index_t& seq_index) const {
                return _x0 == seq_index[0] && _x1 == seq_index[1];
        }
        virtual __inline__ void operator=(const seq_index_t& seq_index) {
                _x0 = seq_index[0];
                _x1 = seq_index[1];
        }
        virtual __inline__ void operator++(int i) {
                _x1++;
        }

protected:
        size_t _x1;
};

class range_t {
public:
        range_t(const index_t& index, size_t length) 
                : index(index), length(length) {
        }

        const index_t& index;
        const size_t length;
};

#endif /* INDEX_HH */
