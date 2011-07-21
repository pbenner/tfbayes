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
        explicit index_t(size_t x0);
        index_t(size_t x0, size_t x1);
        index_t(const index_t& index);

        friend std::ostream& operator<< (std::ostream& o, const index_t& index);

        size_t operator[](size_t i) const;
        size_t operator[](size_t i);

        void operator=(const index_t& index);
        void operator++(int i);
        bool operator<(index_t index) const;
        size_t size() const;

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
