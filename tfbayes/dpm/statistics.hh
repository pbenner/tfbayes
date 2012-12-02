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

#ifndef STATISTICS_HH
#define STATISTICS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

static inline
size_t select_component(size_t k, double log_weights[])
{
        /* log_weights are cumulative */

        const double r = (double)rand()/RAND_MAX;

        if (r == 0.0) {
                return 0;
        }

        const double log_r = log(r) + log_weights[k-1];
        for (size_t i = 0; i < k-1; i++) {
                if (log_r <= log_weights[i]) {
                        return i;
                }
        }
        return k-1;
}

#endif /* STATISTICS_HH */
