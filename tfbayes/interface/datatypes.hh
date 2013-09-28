/* Copyright (C) 2012 Philipp Benner
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

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

// c++ strings
////////////////////////////////////////////////////////////////////////////////

#include <string>

__BEGIN_DECLS

std::string* _cxx_string_alloc(const char* str);
void _cxx_string_free(std::string* ptr);

__END_DECLS

// c++ vector/matrix interface
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cstddef>

#include <tfbayes/utility/linalg.hh>

__BEGIN_DECLS

std::vector<double>* _cxx_vector_alloc(size_t size);
std::matrix<double>* _cxx_matrix_alloc(size_t rows, size_t columns);

double _cxx_vector_read(std::vector<double>* vector, size_t pos);
double _cxx_matrix_read(std::matrix<double>* matrix, size_t row, size_t column);

void _cxx_vector_write(std::vector<double>* vector, size_t pos, double d);
void _cxx_matrix_write(std::matrix<double>* matrix, size_t row, size_t column, double d);

void _cxx_vector_free(std::vector<double>* vector);
void _cxx_matrix_free(std::matrix<double>* matrix);

__END_DECLS

// vector interface
////////////////////////////////////////////////////////////////////////////////

__BEGIN_DECLS

vector_t * _alloc_vector(size_t size);
matrix_t * _alloc_matrix(size_t rows, size_t columns);
void _free_vector(vector_t *v);
void _free_matrix(matrix_t *m);
void _free(void *ptr);

__END_DECLS
