/* Copyright (C) 2010 Philipp Benner
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

#ifndef _FASTLOG_H_
#define _FASTLOG_H_

#include <math.h>

#define LOGTABLE_SIZE 1000000

extern double logtable[LOGTABLE_SIZE];

static inline
double fastlog(size_t a)
{
        return a > LOGTABLE_SIZE ? log(a) : logtable[a];
}

#endif /* _FASTLOG_H_ */
