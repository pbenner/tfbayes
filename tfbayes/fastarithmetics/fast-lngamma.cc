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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cstdlib>

#include <gsl/gsl_sf_gamma.h>

#include <fast-lngamma.hh>
#include <fast-lngamma-table.hh>

double fast_lngamma(double x)
{
        if (x <= 1000.0) {
                // static casting uses a simple trunc() to
                // convert a double to integer, therefore add
                // 0.5 to reduce the numerical error
                size_t n = static_cast<size_t>(x*100.0-1.0+0.5);

                return fast_lngamma_table[n][1];
        }
        else {
                return gsl_sf_lngamma(x);
        }
}
