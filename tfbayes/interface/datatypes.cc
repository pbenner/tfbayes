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

#include <iostream>
#include <tfbayes/interface/datatypes.hh>

using namespace std;

__BEGIN_DECLS

// c++ strings
////////////////////////////////////////////////////////////////////////////////

string* _cxx_string_alloc(const char* str)
{
        string* ptr = new string(str);
        return ptr;
}

void _cxx_string_free(string* ptr)
{
        delete(ptr);
}

const char* _cxx_string_getstr(string* ptr)
{
        return ptr->c_str();
}

// c++ vector/matrix interface
////////////////////////////////////////////////////////////////////////////////

vector<double>* _cxx_vector_alloc(size_t size)
{
        vector<double>* ptr = new vector<double>(size, 0.0);
        return ptr;
}

matrix<double>* _cxx_matrix_alloc(size_t rows, size_t columns)
{
        matrix<double>* ptr = new matrix<double>(rows, columns);
        return ptr;
}

double _cxx_vector_read(vector<double>* vector, size_t pos)
{
        return vector->operator[](pos);
}

double _cxx_matrix_read(matrix<double>* matrix, size_t row, size_t column)
{
        return matrix->operator[](row)[column];
}

void _cxx_vector_write(vector<double>* vector, size_t pos, double d)
{
        vector->operator[](pos) = d;
}

void _cxx_matrix_write(matrix<double>* matrix, size_t row, size_t column, double d)
{
        matrix->operator[](row)[column] = d;
}

void _cxx_vector_free(std::vector<double>* vector)
{
        delete(vector);
}

void _cxx_matrix_free(std::matrix<double>* matrix)
{
        delete(matrix);
}

// vector interface
////////////////////////////////////////////////////////////////////////////////

vector_t * _alloc_vector(size_t size) {
        return alloc_vector(size);
}
matrix_t * _alloc_matrix(size_t rows, size_t columns) {
        return alloc_matrix(rows, columns);
}
void _free_vector(vector_t *v) { free_vector(v); }
void _free_matrix(matrix_t *m) { free_matrix(m); }
void _free(void *ptr)          { free(ptr); }

__END_DECLS
