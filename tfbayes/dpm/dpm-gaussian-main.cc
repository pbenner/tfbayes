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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <gsl/gsl_matrix.h>

#include <tfbayes/interface/datatypes.hh>
#include <tfbayes/dpm/init.hh>
#include <tfbayes/dpm/data-gaussian.hh>
#include <tfbayes/dpm/dpm-gaussian.hh>
#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/utility/linalg.hh>

void sample(
        int samples,
        double alpha,
        matrix_t* _Sigma,
        matrix_t* _Sigma_0,
        vector_t* _mu_0,
        vector_t* _pi)
{
        __dpm_init__();

        gsl_matrix *Sigma    = to_gsl_matrix(_Sigma);
        gsl_matrix *Sigma_0  = to_gsl_matrix(_Sigma_0);
        gsl_vector *mu_0     = to_gsl_vector(_mu_0);
        const size_t cluster = _pi->size;

        data_gaussian_t* _data    = new data_gaussian_t(cluster, (size_t)samples, Sigma, _pi->vec);
        dpm_gaussian_t*  _gdpm    = new dpm_gaussian_t(alpha, Sigma, Sigma_0, mu_0, *_data);
        gibbs_sampler_t* _sampler = new gibbs_sampler_t(*_gdpm, *_gdpm, *_data);

        gsl_matrix_free(Sigma);
        gsl_matrix_free(Sigma_0);
        gsl_vector_free(mu_0);

        _sampler->sample(100, 100);

}

int
main(void) {
        matrix_t* Sigma   = alloc_matrix(2, 2);
        matrix_t* Sigma_0 = alloc_matrix(2, 2);
        vector_t* mu      = alloc_vector(2);
        vector_t* pi      = alloc_vector(11);

        Sigma->mat[0][0] = 0.01;
        Sigma->mat[0][1] = 0.005;
        Sigma->mat[1][0] = 0.005;
        Sigma->mat[1][1] = 0.01;

        Sigma_0->mat[0][0] = 10.0;
        Sigma_0->mat[0][1] = 0.1;
        Sigma_0->mat[1][0] = 0.1;
        Sigma_0->mat[1][1] = 10.0;

        mu->vec[0] = 0;
        mu->vec[1] = 0;

        pi->vec[0]  = 0.27693787;
        pi->vec[1]  = 0.06001137;
        pi->vec[2]  = 0.10600994;
        pi->vec[3]  = 0.00997665;
        pi->vec[4]  = 0.02111005;
        pi->vec[5]  = 0.0120215;
        pi->vec[6]  = 0.04216835;
        pi->vec[7]  = 0.06136474;
        pi->vec[8]  = 0.05276006;
        pi->vec[9]  = 0.29406385;
        pi->vec[10] = 0.06357562;

        sample(1000, 1, Sigma, Sigma_0, mu, pi);
}
