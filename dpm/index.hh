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

#include <iostream>

// index_t and range_t
////////////////////////////////////////////////////////////////////////////////

class index_t {
public:
        __inline__ explicit index_t(size_t x0) : _x0(x0), _size(1) {
                _x0 = x0;
        }
        __inline__ index_t(size_t x0, size_t x1) : _x0(x0), _x1(x1), _size(2) {
        }
        __inline__ index_t(const index_t& index) : _x0(index._x0), _x1(index._x1), _size(index._size) {
        }

        friend std::ostream& operator<< (std::ostream& o, const index_t& index);

        __inline__ const size_t& operator[](size_t i) const {
                return i == 0 ? _x0 : _x1;
        }
        __inline__ size_t& operator[](size_t i) {
                return i == 0 ? _x0 : _x1;
        }

        __inline__ void operator=(const index_t& index) {
                _x0 = index[0];
                _x1 = index[1];
        }
        __inline__ void operator++(int i) {
                _size == 1 ? _x0++ : _x1++;
        }
        __inline__ bool operator<(index_t index) const {
                if (_size == 1) {
                        return _x0 < index[0];
                }
                else {
                        return _x0 < index[0] || (_x0 == index[0] && _x1 < index[1]);
                }
        }
        __inline__ size_t size() const {
                return _size;
        }


private:
        size_t _x0;
        size_t _x1;
        const size_t _size;
};

class range_t {
public:
        range_t(const index_t& from, const index_t& to) 
                : from(from), to(to){
        }

        const index_t from;
        const index_t to;
};

#endif /* INDEX_HH */
