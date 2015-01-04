/* Copyright (C) 2010-2013 Philipp Benner
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

#ifndef __TFBAYES_UTILITY_LINALG_HH__
#define __TFBAYES_UTILITY_LINALG_HH__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <numeric> /* accumulate */
#include <vector>
#include <cstddef>

/* since this is the only position where c++0x might be useful
 * in the code we try to stick with boost here */
//#include <type_traits>
#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/enable_if.hpp>

#include <gsl/gsl_matrix.h>

// c++ matrix
////////////////////////////////////////////////////////////////////////////////

namespace std {
        template <typename T>
        class matrix : public vector<vector<T> > {
        public:
                matrix()
                        : vector<vector<T> >()
                        {}
                matrix(size_t rows, size_t columns)
                        : vector<vector<T> >(rows, vector<T>(columns, 0.0))
                        {}
                matrix(size_t rows, size_t columns, T init)
                        : vector<vector<T> >(rows, vector<T>(columns, init))
                        {}
                template <typename InputIterator>
                matrix (InputIterator first,
                        InputIterator last,
                        const allocator<vector<double> >& a = allocator<vector<double> >(),
                        typename boost::disable_if<boost::is_integral<InputIterator> >::type* = 0)
                        : vector<vector<T> >(first, last, a) { }
                matrix<T> transpose(T init = 0.0) {
                        matrix<T> m(columns(), rows(), init);
                        for (size_t i = 0; i < columns(); i++) {
                                for (size_t j = 0; j < rows(); j++) {
                                        m[i][j] = (*this)[j][i]; 
                                }
                        }
                        return m;
                }
                size_t rows() const {
                        return matrix<T>::size();
                }
                size_t columns() const {
                        return rows() == 0 ? 0 : operator[](0).size();
                }
                using vector<vector<T> >::operator[];
        };
}

#include <boost/foreach.hpp>

template <typename T>
std::vector<T> normalize(const std::vector<T>& container)
{
        std::vector<T> tmp(container);
        T sum = std::accumulate(tmp.begin(), tmp.end(), static_cast<T>(0));

        BOOST_FOREACH(T &n, tmp) {
                n /= sum;
        }
        return tmp;
}

// gsl vector/matrix conversion
////////////////////////////////////////////////////////////////////////////////

static inline
gsl_vector * to_gsl_vector(const std::vector<double>& vector)
{
        gsl_vector *v = gsl_vector_alloc(vector.size());

        for(size_t i = 0; i < vector.size(); i++) {
                gsl_vector_set(v, i, vector[i]);
        }
        return v;
}

static inline
std::vector<double> from_gsl_vector(const gsl_vector * vector)
{
        std::vector<double> v(vector->size, 0);

        for(size_t i = 0; i < vector->size; i++) {
                v[i] = gsl_vector_get(vector, i);
        }
        return v;
}

static inline
gsl_matrix * to_gsl_matrix(const std::matrix<double>& matrix)
{
        gsl_matrix *m = gsl_matrix_alloc(matrix.rows(), matrix.columns());

        for(size_t i = 0; i < matrix.rows(); i++) {
                for(size_t j = 0; j < matrix.columns(); j++) {
                        gsl_matrix_set(m, i, j, matrix[i][j]);
                }
        }
        return m;
}

static inline
std::matrix<double> from_gsl_matrix(const gsl_matrix * matrix)
{
        std::matrix<double> m(matrix->size1, matrix->size2);

        for(size_t i = 0; i < matrix->size1; i++) {
                for(size_t j = 0; j < matrix->size2; j++) {
                        m[i][j] = gsl_matrix_get(matrix, i, j);
                }
        }
        return m;
}

#endif /* __TFBAYES_UTILITY_LINALG_HH__ */
