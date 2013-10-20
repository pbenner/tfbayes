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

#ifndef _LOGARITHMETIC_H_
#define _LOGARITHMETIC_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cstdlib>
#include <cmath>
#include <limits>

/* Log Sum of Exponentials Algorithm */

static inline
double logadd(double a, double b)
{
        if (a < b) return a == -std::numeric_limits<double>::infinity() ? b : b + log1p(exp(a-b));
        else       return b == -std::numeric_limits<double>::infinity() ? a : a + log1p(exp(b-a));
}

static inline
double logsub(double a, double b)
{
        return b == -std::numeric_limits<double>::infinity() ? a : a + log(1-exp(b-a));
}

#endif /* _LOGARITHMETIC_H_ */
