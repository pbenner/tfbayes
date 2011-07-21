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

#include <index.hh>

using namespace std;

index_t::index_t(size_t x0) : _x0(x0), _size(1) {
        _x0 = x0;
}

index_t::index_t(size_t x0, size_t x1) : _x0(x0), _x1(x1), _size(2) {
}

index_t::index_t(const index_t& index) : _x0(index._x0), _x1(index._x1), _size(index._size) {
}

size_t
index_t::operator[](size_t i) const {
        return i == 0 ? _x0 : _x1;
}

size_t
index_t::operator[](size_t i) {
        return i == 0 ? _x0 : _x1;
}

void
index_t::operator=(const index_t& index) {
        _x0 = index[0];
        _x1 = index[1];
}

void
index_t::operator++(int i) {
        _size == 1 ? _x0++ : _x1++;
}

bool
index_t::operator<(index_t index) const
{
        if (_size == 1) {
                return _x0 < index[0];
        }
        else {
                return _x0 < index[0] || (_x0 == index[0] && _x1 < index[1]);
        }
}

size_t
index_t::size() const {
        return _size;
}
