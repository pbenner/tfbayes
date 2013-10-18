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

#ifndef DATA_HH
#define DATA_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include <tfbayes/dpm/index.hh>
#include <tfbayes/dpm/datatypes.hh>
#include <tfbayes/utility/clonable.hh>

// data_t and iterator_t
////////////////////////////////////////////////////////////////////////////////

template <typename T> class data_i;
template <typename T> class data_t;
template <typename T> class sequence_data_t;

template <typename T>
std::ostream&
operator<< (std::ostream& o, const data_t<T>& sd) {
        for (size_t i = 0; i < sd.size(); i++) {
                o << sd[i] << " ";
        }
        o << std::endl;

        return o;
}

template <typename T>
std::ostream&
operator<< (std::ostream& o, const sequence_data_t<T>& sd) {
        for (size_t i = 0; i < sd.size(); i++) {
                for (size_t j = 0; j < sd[i].size(); j++) {
                        o << sd[i][j] << " ";
                }
                o << std::endl;
        }
        return o;
}

template <typename T>
class iterator_t : public std::iterator<std::forward_iterator_tag, T>
{
public:
        iterator_t(data_i<T>& data, const index_i& index, size_t length)
                : _data(data), _index(index), _pos(*index.clone()), _length(length), _i(0) {
        }
        ~iterator_t() {
                delete(&_pos);
        }

        const T& operator*() const {
                return _data[_pos];
        }
        T& operator*() {
                return _data[_pos];
        }
        bool operator++(int i) {
                if (_i+1 < _length) {
                        _pos++; _i++;
                        return true;
                }
                return false;
        }
        bool operator==(const iterator_t& it) const {
                return (_index == it._index) &&
                       (_pos   == it._pos  );
        }
        bool operator!=(const iterator_t& it) const {
                return (_index != it._index) ||
                       (_pos   != it._pos  );
        }
protected:
        data_i<T>& _data;
        const index_i& _index;
        index_i& _pos;
        size_t _length, _i;
};

template <typename T>
class const_iterator_t : public std::iterator<std::forward_iterator_tag, const T>
{
public:
        const_iterator_t(const data_i<T>& data, const index_i& index, size_t length)
                : _data(data), _index(index), _pos(*index.clone()), _length(length), _i(0) {
        }
        ~const_iterator_t() {
                delete(&_pos);
        }

        const T& operator*() const {
                return _data[_pos];
        }
        bool operator++(int i) {
                if (_i+1 < _length) {
                        _pos++; _i++;
                        return true;
                }
                return false;
        }
        bool operator==(const const_iterator_t& it) const {
                return (_index == it._index) &&
                       (_pos   == it._pos  );
        }
        bool operator!=(const const_iterator_t& it) const {
                return (_index != it._index) ||
                       (_pos   != it._pos  );
        }
protected:
        const data_i<T>& _data;
        const index_i& _index;
        index_i& _pos;
        size_t _length, _i;
};

template <typename T>
class data_i : public clonable {
public:
        virtual data_i<T>* clone() const = 0;
        virtual data_i<T>& operator=(const data_i<T>& data) {
                std::cerr << "Error: data_i<T> assignment operator called"
                          << std::endl;
                exit(EXIT_FAILURE);
        }
        virtual const_iterator_t<T> operator[](const range_t& range) const = 0;
        virtual iterator_t<T> operator[](const range_t& range) = 0;
        virtual const T& operator[](const index_i& index) const = 0;
        virtual T& operator[](const index_i& index) = 0;
};

template <typename T>
class data_t : public data_i<T>, public std::vector<T>
{
public:
        using std::vector<T>::size;
        using std::vector<T>::push_back;
        using std::vector<T>::operator[];
        using std::vector<T>::operator=;

        data_t() : std::vector<T>() {
        }
        data_t(size_t n, T init)
                : std::vector<T>(n, init) {
        }
        data_t(const std::vector<T>& data)
                : std::vector<T>(data) {
        }
        virtual ~data_t() {}

        friend void swap(data_t<T>& first, data_t<T>& second) {
                using std::swap;
                swap(static_cast<std::vector<T>&>(first),
                     static_cast<std::vector<T>&>(second));
        }
        virtual data_t<T>* clone() const {
                return new data_t<T>(*this);
        }
        virtual data_t<T>& operator=(const data_i<T>& data) {
                std::vector<T>::operator=(
                        static_cast<const data_t<T>&>(data));
                return *this;
        }
        virtual const_iterator_t<T> operator[](const range_t& range) const {
                return const_iterator_t<T>(*this, range.index(), range.length());
        }
        virtual iterator_t<T> operator[](const range_t& range) {
                return iterator_t<T>(*this, range.index(), range.length());
        }
        virtual const T& operator[](const index_i& index) const {
                return std::vector<T>::operator[](index[0]);
        }
        virtual T& operator[](const index_i& index) {
                return std::vector<T>::operator[](index[0]);
        }
        friend std::ostream& operator<< <> (std::ostream& o, const data_t<T>& sd);
};

template <typename T>
class sequence_data_t : public data_i<T>, public std::vector<std::vector<T> >
{
public:
        using std::vector<std::vector<T> >::size;
        using std::vector<std::vector<T> >::push_back;
        using std::vector<std::vector<T> >::operator[];
        using std::vector<std::vector<T> >::operator=;

        sequence_data_t() : std::vector<std::vector<T> >() {
        }
        sequence_data_t(const std::vector<size_t> n, const T init)
                : std::vector<std::vector<T> >() {
                for (size_t i = 0; i < n.size(); i++) {
                        push_back(std::vector<T>(n[i], init));
                }
        }
        sequence_data_t(const std::vector<std::vector<T> >& data)
                : std::vector<std::vector<T> >(data) {
        }
        virtual ~sequence_data_t() {
        }

        void swap(sequence_data_t<T>& first, sequence_data_t<T>& second) {
                using std::swap;
                swap(static_cast<std::vector<std::vector<T> >&>(first),
                     static_cast<std::vector<std::vector<T> >&>(second));
        }
        virtual sequence_data_t<T>* clone() const {
                return new sequence_data_t<T>(*this);
        }
        virtual sequence_data_t<T>& operator=(const data_i<T>& data) {
                std::vector<std::vector<T> >::operator=(
                        static_cast<const sequence_data_t<T>&>(data));
                return *this;
        }
        virtual const_iterator_t<T> operator[](const range_t& range) const {
                return const_iterator_t<T>(*this, range.index(), range.length());
        }
        virtual iterator_t<T> operator[](const range_t& range) {
                return iterator_t<T>(*this, range.index(), range.length());
        }
        virtual const T& operator[](const index_i& index) const {
                return std::vector<std::vector<T> >::operator[](index[0])[index[1]];
        }
        virtual T& operator[](const index_i& index) {
                return std::vector<std::vector<T> >::operator[](index[0])[index[1]];
        }
        virtual size_t size(size_t i) const {
                return operator[](i).size();
        }
        virtual std::vector<size_t> sizes() const {
                std::vector<size_t> lengths;
                for (size_t i = 0; i < size(); i++) {
                        lengths.push_back(operator[](i).size());
                }
                return lengths;
        }
        friend std::ostream& operator<< <> (std::ostream& o, const sequence_data_t<T>& sd);
};

#endif /* DATA_HH */
