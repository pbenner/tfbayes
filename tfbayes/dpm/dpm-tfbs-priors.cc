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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <dpm-tfbs.hh>

double
DpmTfbs::py_prior(Cluster& cluster)
{
        if (cluster.size() == 0) {
                return log(alpha + discount*(mixture_components()-1)) - log(num_tfbs + alpha);
        }
        else {
                return log(cluster.size()-discount) - log(num_tfbs + alpha);
        }
}

double
DpmTfbs::uniform_prior(Cluster& cluster)
{
        double K = mixture_components()-1;

        if (cluster.size() == 0) {
                return log(alpha) - log(alpha + K);
        }
        else {
                return -log(alpha + K);
        }
}

double
DpmTfbs::poppe_prior(Cluster& cluster)
{
        double K = mixture_components()-1;
        double N = num_tfbs;

        if (K == 0 && cluster.size() == 0) {
                return 0;
        }
        if (cluster.size() == 0) {
                if (K == 1.0) {
                        return -log(N+1.0);
                }
                else {
                        return log(K*(K-1.0)/(N*(N+1.0)));
                }
        }
        else {
                if (K == 1.0) {
                        return log(N/(N+1.0));
                }
                else {
                        return log((cluster.size()+1.0)/(N+1.0) * (N-K+1.0)/N);
                }
        }
}
