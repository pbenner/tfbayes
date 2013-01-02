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
#include <tfbayes/config.h>
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

static inline
size_t select_max_component(size_t k, double log_weights[])
{
        /* log_weights are cumulative, so we need to normalize by log_weights[k-1] */

        /* resulting cluster number */
        size_t result = 0;
        /* difference between two successive log weights */
        double diff   = 0.0;
        /* we need to store the last value because zero is not
         * represented in the log weights array */
        double last   = 0.0;

        for (size_t i = 0; i < k; i++) {
                /* if the increase for this cluster is higher than any
                 * before */
                if (diff < exp(log_weights[i] - log_weights[k-1]) - last) {
                        /* store the new increase in probability */
                        diff   = exp(log_weights[i] - log_weights[k-1]) - last;
                        /* and the cluster number */
                        result = i;
                }
                /* also save the value of the last weight */
                last = exp(log_weights[i] - log_weights[k-1]);
        }
        return result;
}

#endif /* STATISTICS_HH */
