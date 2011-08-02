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

#include <iostream>

#include <interface.hh>
#include <dpm-gaussian-interface.hh>

int
main(void) {
        Bayes::Matrix* Sigma   = Bayes::allocMatrix(2, 2);
        Bayes::Matrix* Sigma_0 = Bayes::allocMatrix(2, 2);
        Bayes::Vector* mu      = Bayes::allocVector(2);
        Bayes::Vector* pi      = Bayes::allocVector(11);

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

        _dpm_gaussian_init(1000, 1, Sigma, Sigma_0, mu, pi);
        _dpm_gaussian_sample(100,100);
}
