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

#ifndef __TFBAYES_DPM_DATA_HH__
#define __TFBAYES_DPM_DATA_HH__

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

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

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
class data_i : public virtual clonable {
public:
        class iterator_t : public std::iterator<std::forward_iterator_tag, T>
        {
        public:
                iterator_t(data_i<T>& data, const index_t& index, size_t length)
                        : m_data  (&data)
                        , m_index (index)
                        , m_pos   (index)
                        , m_length(length)
                        , m_i     (0) {
                }

                const T& operator*() const {
                        return (*m_data)[m_pos];
                }
                T& operator*() {
                        return (*m_data)[m_pos];
                }
                bool operator++(int i) {
                        if (m_i+1 < m_length) {
                                m_pos++; m_i++;
                                return true;
                        }
                        return false;
                }
                bool operator==(const iterator_t& it) const {
                        return (m_index == it.m_index) &&
                                (m_pos   == it.m_pos  );
                }
                bool operator!=(const iterator_t& it) const {
                        return (m_index != it.m_index) ||
                                (m_pos   != it.m_pos  );
                }
        protected:
                data_i<T>* m_data;
                index_t m_index;
                index_t m_pos;
                size_t m_length, m_i;
        };

        class const_iterator_t : public std::iterator<std::forward_iterator_tag, const T>
        {
        public:
                const_iterator_t(const data_i<T>& data, const index_t& index, size_t length)
                        : m_data  (&data)
                        , m_index (index)
                        , m_pos   (index)
                        , m_length(length)
                        , m_i     (0) {
                }

                const T& operator*() const {
                        return (*m_data)[m_pos];
                }
                bool operator++(int i) {
                        if (m_i+1 < m_length) {
                                m_pos++; m_i++;
                                return true;
                        }
                        return false;
                }
                bool operator==(const const_iterator_t& it) const {
                        return (m_index == it.m_index) &&
                                (m_pos   == it.m_pos  );
                }
                bool operator!=(const const_iterator_t& it) const {
                        return (m_index != it.m_index) ||
                                (m_pos   != it.m_pos  );
                }
        protected:
                const data_i<T>* m_data;
                index_t m_index;
                index_t m_pos;
                size_t m_length, m_i;
        };

        virtual data_i<T>* clone() const = 0;
        virtual const_iterator_t operator[](const range_t& range) const = 0;
        virtual iterator_t operator[](const range_t& range) = 0;
        virtual const T& operator[](const index_t& index) const = 0;
        virtual T& operator[](const index_t& index) = 0;
};

template <typename T>
class data_t : public data_i<T>, public std::vector<T>
{
public:
        typedef std::vector<T> base_t;
        typedef typename data_i<T>::iterator_t iterator_t;
        typedef typename data_i<T>::const_iterator_t const_iterator_t;
        using base_t::size;
        using base_t::push_back;
        using base_t::operator[];

        data_t() : base_t() {
        }
        data_t(size_t n, T init)
                : base_t(n, init) {
        }
        data_t(const data_t& data)
                : base_t(data)
                { }

        virtual data_t<T>* clone() const {
                return new data_t<T>(*this);
        }

        data_t& operator=(const data_t& data) {
                base_t::operator=(data);
                return *this;
        }
        virtual inline const_iterator_t operator[](const range_t& range) const GCC_ATTRIBUTE_HOT {
                return const_iterator_t(*this, range.index(), range.length());
        }
        virtual inline iterator_t operator[](const range_t& range) GCC_ATTRIBUTE_HOT {
                return iterator_t(*this, range.index(), range.length());
        }
        virtual inline const T& operator[](const index_t& index) const GCC_ATTRIBUTE_HOT {
                return base_t::operator[](index[0]);
        }
        virtual inline T& operator[](const index_t& index) GCC_ATTRIBUTE_HOT {
                return base_t::operator[](index[0]);
        }
        friend std::ostream& operator<< <> (std::ostream& o, const data_t<T>& sd);
private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) {
                ar & boost::serialization::base_object<base_t >(*this);
        }
};

template <typename T>
class sequence_data_t : public data_i<T>, public std::vector<std::vector<T> >
{
public:
        typedef std::vector<std::vector<T> > base_t;
        typedef typename data_i<T>::iterator_t iterator_t;
        typedef typename data_i<T>::const_iterator_t const_iterator_t;
        using base_t::size;
        using base_t::push_back;
        using base_t::operator[];
        using base_t::operator=;

        sequence_data_t() : base_t() {
        }
        sequence_data_t(const std::vector<size_t> n, const T init)
                : base_t() {
                for (size_t i = 0; i < n.size(); i++) {
                        push_back(std::vector<T>(n[i], init));
                }
        }
        sequence_data_t(const sequence_data_t& sequence_data)
                : base_t(sequence_data)
                { }

        virtual sequence_data_t<T>* clone() const {
                return new sequence_data_t<T>(*this);
        }

        sequence_data_t& operator=(const sequence_data_t& sequence_data) {
                base_t::operator=(sequence_data);
                return *this;
        }
        virtual inline const_iterator_t operator[](const range_t& range) const GCC_ATTRIBUTE_HOT {
                return const_iterator_t(*this, range.index(), range.length());
        }
        virtual inline iterator_t operator[](const range_t& range) GCC_ATTRIBUTE_HOT {
                return iterator_t(*this, range.index(), range.length());
        }
        virtual inline const T& operator[](const index_t& index) const GCC_ATTRIBUTE_HOT {
                return base_t::operator[](index[0])[index[1]];
        }
        virtual inline T& operator[](const index_t& index) GCC_ATTRIBUTE_HOT {
                return base_t::operator[](index[0])[index[1]];
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
private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) {
                ar & boost::serialization::base_object<base_t>(*this);
        }
};

#endif /* __TFBAYES_DPM_DATA_HH__ */
