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

#include <interface.hh>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

Simple::Vector * _allocVector(int size)              { return Simple::allocVector(size); }
void            _freeVector(Simple::Vector *v)       { Simple::freeVector(v); }
Simple::Matrix * _allocMatrix(int rows, int columns) { return Simple::allocMatrix(rows, columns); }
void            _freeMatrix(Simple::Matrix *m)       { Simple::freeMatrix(m); }
void            _free(void *ptr)                     { free(ptr); }

__END_DECLS
