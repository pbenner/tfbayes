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

#ifndef DATATYPES_HH
#define DATATYPES_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <vector>
#include <string>

// index_t and range_t
////////////////////////////////////////////////////////////////////////////////

class index_t {
public:
        explicit index_t(size_t x0) : _x0(x0), _size(1) {
                _x0 = x0;
        }
        index_t(size_t x0, size_t x1) : _x0(x0), _x1(x1), _size(2) {
        }
        index_t(const index_t& index) : _x0(index._x0), _x1(index._x1), _size(index._size) {
        }

        friend std::ostream& operator<< (std::ostream& o, const index_t& index);

        inline size_t operator[](size_t i) const { return i == 0 ? _x0 : _x1; }
        inline size_t operator[](size_t i)       { return i == 0 ? _x0 : _x1; }

        void operator=(const index_t& index) {
                _x0 = index[0];
                _x1 = index[1];
        }

        void operator++(int i) {
                _size == 1 ? _x0++ : _x1++;
        }

        bool operator<(index_t index) const {
                if (_size == 1) {
                        return _x0 < index[0];
                }
                else {
                        return _x0 < index[0] || (_x0 == index[0] && _x1 < index[1]);
                }
        }

        size_t size() const {
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
        iterator_t(data_t<T>& data, const index_t& from, const index_t& to)
                : _data(data), _from(from), _to(to), _pos(from) {
        }

        const T& operator*() const {
                return _data[_pos];
        }
              T& operator*() {
                return _data[_pos];
        }

        void reset() {
                _pos = _from;
        }

        const bool operator++(int i) {
                if (_pos < _to) {
                        _pos++;
                        return true;
                }
                return false;
        }

private:
              data_t<T>& _data;
        const index_t& _from;
        const index_t& _to;
        index_t _pos;
};

template <typename T>
class const_iterator_t
{
public:
        const_iterator_t(const data_t<T>& data, const index_t& from, const index_t& to)
                : _data(data), _from(from), _to(to), _pos(from) {
        }

        const T& operator*() const {
                return _data[_pos];
        }

        void reset() {
                _pos = _from;
        }

        const bool operator++(int i) {
                if (_pos < _to) {
                        _pos++;
                        return true;
                }
                return false;
        }

private:
        const data_t<T>& _data;
        const index_t& _from;
        const index_t& _to;
        index_t _pos;
};

template <typename T>
class data_t
{
public:
        data_t() : _data() {
        }
        data_t(size_t n, T init) : _data(n, init) {
        }
        data_t(const std::vector<T>& data) : _data(data) {
        }

        friend std::ostream& operator<< <> (std::ostream& o, const data_t<T>& sd);

        virtual const_iterator_t<T> operator[](const range_t& range) const {
                return const_iterator_t<T>(*this, range.from, range.to);
        }
        virtual iterator_t<T> operator[](const range_t& range) {
                return iterator_t<T>(*this, range.from, range.to);
        }
        virtual const T& operator[](const index_t& index) const {
                return _data[index[0]];
        }
        virtual T& operator[](const index_t& index) {
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

        friend std::ostream& operator<< <> (std::ostream& o, const sequence_data_t<T>& sd);

        virtual const_iterator_t<T> operator[](const range_t& range) const {
                return const_iterator_t<T>((const data_t<T>&)*this, range.from, range.to);
        }
        virtual iterator_t<T> operator[](const range_t& range) {
                return iterator_t<T>((data_t<T>&)*this, range.from, range.to);
        }
        virtual const T& operator[](const index_t& index) const {
                return _data[index[0]][index[1]];
        }

        virtual T& operator[](const index_t& index) {
                return _data[index[0]][index[1]];
        }

protected:
        std::vector<std::vector<T> > _data;
};

// cluster structures
////////////////////////////////////////////////////////////////////////////////

typedef ssize_t cluster_tag_t;

typedef enum {
        cluster_event_empty, cluster_event_nonempty,
        cluster_event_add_word, cluster_event_remove_word
} cluster_event_t;

typedef struct {
        std::vector<std::vector<double> > switches;
        std::vector<std::vector<double> > likelihood;
        std::vector<std::vector<size_t> > components;
} sampling_history_t;

typedef std::vector<std::vector<double> > posterior_t;

#endif /* DATATYPES_HH */
