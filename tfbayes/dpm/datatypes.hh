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

#include <clonable.hh>
#include <index.hh>
#include <data.hh>
#include <graph.hh>

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
        virtual ~sequence_data_t() {}

        virtual sequence_data_t* clone() const {
                return new sequence_data_t(*this);
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
typedef ssize_t model_tag_t;

typedef enum {
        cluster_event_empty, cluster_event_nonempty,
        cluster_event_add_word, cluster_event_remove_word
} cluster_event_t;

typedef struct {
        std::vector<std::vector<double> > switches;
        std::vector<std::vector<double> > likelihood;
        std::vector<std::vector<size_t> > components;
} sampling_history_t;

typedef struct {
        std::vector<std::vector<double> > probabilities;
        std::vector<std::string> hypergraph;
        Graph graph;
} posterior_t;

#endif /* DATATYPES_HH */
