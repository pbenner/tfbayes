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

template <typename T> class data_t;
template <typename T> class sequence_data_t;

template <typename T>
std::ostream&
operator<< (std::ostream& o, const data_t<T>& sd) {
        for (size_t i = 0; i < sd._data.size(); i++) {
                o << sd._data[i] << " ";
        }
        o << std::endl;

        return o;
}

template <typename T>
std::ostream&
operator<< (std::ostream& o, const sequence_data_t<T>& sd) {
        for (size_t i = 0; i < sd._data.size(); i++) {
                for (size_t j = 0; j < sd._data[i].size(); j++) {
                        o << sd._data[i][j] << " ";
                }
                o << std::endl;
        }
        return o;
}

template <typename T>
class iterator_t
{
public:
        iterator_t(data_t<T>& data, const index_i& index, size_t length)
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

        void reset() {
                _pos = _index;
                _i   = 0;
        }

        bool operator++(int i) {
                if (_i+1 < _length) {
                        _pos++; _i++;
                        return true;
                }
                return false;
        }

private:
        data_t<T>& _data;
        const index_i& _index;
        index_i& _pos;
        size_t _length, _i;
};

template <typename T>
class const_iterator_t
{
public:
        const_iterator_t(const data_t<T>& data, const index_i& index, size_t length)
                : _data(data), _index(index), _pos(*index.clone()), _length(length), _i(0) {
        }
        ~const_iterator_t() {
                delete(&_pos);
        }

        const T& operator*() const {
                return _data[_pos];
        }

        void reset() {
                _pos = _index;
                _i   = 0;
        }

        bool operator++(int i) {
                if (_i+1 < _length) {
                        _pos++; _i++;
                        return true;
                }
                return false;
        }

private:
        const data_t<T>& _data;
        const index_i& _index;
        index_i& _pos;
        size_t _length, _i;
};

template <typename T>
class data_t : public clonable
{
public:
        data_t() : _data() {
        }
        data_t(size_t n, T init) : _data(n, init) {
        }
        data_t(const std::vector<T>& data) : _data(data) {
        }
        virtual ~data_t() {}

        virtual data_t* clone() const {
                return new data_t(*this);
        }

        friend void swap(data_t<T>& first, data_t<T>& second) {
                std::swap(first._data, second._data);
        }

        friend std::ostream& operator<< <> (std::ostream& o, const data_t<T>& sd);

        virtual data_t& operator=(const data_t& data) {
                data_t<T> tmp(data);
                swap(*this, tmp);
                return *this;
        }
        virtual const_iterator_t<T> operator[](const range_t& range) const {
                return const_iterator_t<T>(*this, range.index(), range.length());
        }
        virtual iterator_t<T> operator[](const range_t& range) {
                return iterator_t<T>(*this, range.index(), range.length());
        }
        virtual const T& operator[](const index_i& index) const {
                return _data[index[0]];
        }
        virtual T& operator[](const index_i& index) {
                return _data[index[0]];
        }

protected:
        std::vector<T> _data;
};

template <typename T>
class sequence_data_t : public data_t<T>
{
public:
        sequence_data_t() : _data() {
        }
        sequence_data_t(const std::vector<size_t> n, const T init) : _data() {
                for (size_t i = 0; i < n.size(); i++) {
                        _data.push_back(std::vector<T>(n[i], init));
                }
        }
        sequence_data_t(const std::vector<std::vector<T> >& data) : _data(data) {
        }
        virtual ~sequence_data_t() {
        }

        virtual sequence_data_t* clone() const {
                return new sequence_data_t(*this);
        }

        friend void swap(sequence_data_t<T>& first, sequence_data_t<T>& second) {
                std::swap(first._data, second._data);
        }

        friend std::ostream& operator<< <> (std::ostream& o, const sequence_data_t<T>& sd);

        virtual void push_back(const std::vector<T>& sequence) {
                _data.push_back(sequence);
        }
        virtual sequence_data_t& operator=(const data_t<T>& data) {
                sequence_data_t<T> tmp(static_cast<const sequence_data_t<T>&>(data));
                swap(*this, tmp);
                return *this;
        }
        virtual const_iterator_t<T> operator[](const range_t& range) const {
                return const_iterator_t<T>(*this, range.index(), range.length());
        }
        virtual iterator_t<T> operator[](const range_t& range) {
                return iterator_t<T>(*this, range.index(), range.length());
        }
        virtual const T& operator[](const index_i& index) const {
                return _data[index[0]][index[1]];
        }
        virtual T& operator[](const index_i& index) {
                return _data[index[0]][index[1]];
        }
        virtual const std::vector<T>& operator[](size_t i) const {
                return _data[i];
        }
        virtual std::vector<T>& operator[](size_t i) {
                return _data[i];
        }
        virtual size_t size() const {
                return _data.size();
        }
        virtual size_t size(size_t i) const {
                return _data[i].size();
        }
        virtual std::vector<size_t> sizes() const {
                std::vector<size_t> lengths;
                for (size_t i = 0; i < size(); i++) {
                        lengths.push_back(_data[i].size());
                }
                return lengths;
        }

protected:
        std::vector<std::vector<T> > _data;
};

#endif /* DATA_HH */
