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

__BEGIN_C_REGION;

Bayes::Vector * _allocVector(int size)              { return Bayes::allocVector(size); }
void            _freeVector(Bayes::Vector *v)       { Bayes::freeVector(v); }
Bayes::Matrix * _allocMatrix(int rows, int columns) { return Bayes::allocMatrix(rows, columns); }
void            _freeMatrix(Bayes::Matrix *m)       { Bayes::freeMatrix(m); }
void            _free(void *ptr)                    { free(ptr); }

__END_C_REGION;
