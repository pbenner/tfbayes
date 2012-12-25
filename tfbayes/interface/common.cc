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

#include <tfbayes/interface/common.hh>

__BEGIN_DECLS

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
