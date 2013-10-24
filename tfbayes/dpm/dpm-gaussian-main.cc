/* Copyright (C) 2011-2013 Philipp Benner
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

#include <iostream>

#include <gsl/gsl_matrix.h>

#include <tfbayes/dpm/data-gaussian.hh>
#include <tfbayes/dpm/dpm-gaussian.hh>
#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/utility/linalg.hh>

using namespace std;

void sample(
        size_t samples,
        double alpha,
        const matrix<double> Sigma,
        const matrix<double> Sigma_0,
        const vector<double> mu_0,
        const vector<double> pi)
{
        data_gaussian_t data(samples, Sigma, pi);
        dpm_gaussian_t  gdpm(alpha, Sigma, Sigma_0, mu_0, data);
        gibbs_sampler_t sampler(gdpm, data);

        sampler(100, 100);
}

int
main(void) {
        matrix<double> Sigma   (2,2,0);
        matrix<double> Sigma_0 (2,2,0);
        vector<double> mu_0    ( 2,0);
        vector<double> pi      (11,0);

        Sigma[0][0] = 0.01;
        Sigma[0][1] = 0.005;
        Sigma[1][0] = 0.005;
        Sigma[1][1] = 0.01;

        Sigma_0[0][0] = 10.0;
        Sigma_0[0][1] = 0.1;
        Sigma_0[1][0] = 0.1;
        Sigma_0[1][1] = 10.0;

        mu_0[0] = 0;
        mu_0[1] = 0;

        pi[0]  = 0.27693787;
        pi[1]  = 0.06001137;
        pi[2]  = 0.10600994;
        pi[3]  = 0.00997665;
        pi[4]  = 0.02111005;
        pi[5]  = 0.0120215;
        pi[6]  = 0.04216835;
        pi[7]  = 0.06136474;
        pi[8]  = 0.05276006;
        pi[9]  = 0.29406385;
        pi[10] = 0.06357562;

        sample(1000, 1, Sigma, Sigma_0, mu_0, pi);
}
