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

#include <vector>

#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include "statistics.hh"

using namespace std;

gsl_rng* _r;

BivariateNormal::BivariateNormal(
        const gsl_matrix* cov,
        const gsl_vector* mu)
{
        update(cov, mu);
}

void BivariateNormal::update(
        const gsl_matrix* cov,
        const gsl_vector* mu)
{
        this->mu_x    = gsl_vector_get(mu, 0);
        this->mu_y    = gsl_vector_get(mu, 1);
        this->sigma_x = sqrt(gsl_matrix_get(cov, 0, 0));
        this->sigma_y = sqrt(gsl_matrix_get(cov, 1, 1));
        this->rho     = gsl_matrix_get(cov, 0, 1)/(sigma_x*sigma_y);
}

BivariateNormal::BivariateNormal(
        const vector<double>& mu,
        const vector<double>& sigma,
        double rho)
{
        update(mu, sigma, rho);
}

void BivariateNormal::update(
        const vector<double>& mu,
        const vector<double>& sigma,
        double rho)
{
        this->mu_x    = mu[0];
        this->mu_y    = mu[1];
        this->sigma_x = sigma[0];
        this->sigma_y = sigma[1];
        this->rho     = rho;
}

BivariateNormal::BivariateNormal(
        double mu_x, double mu_y,
        double sigma_x, double sigma_y,
        double rho)
{
        update(mu_x, mu_y, sigma_x, sigma_y, rho);
}

void BivariateNormal::update(
        double mu_x, double mu_y,
        double sigma_x, double sigma_y,
        double rho)
{
        this->mu_x    = mu_x;
        this->mu_y    = mu_y;
        this->sigma_x = sigma_x;
        this->sigma_y = sigma_y;
        this->rho     = rho;
}

BivariateNormal::BivariateNormal()
{
        update(0, 0, 0, 0, 0);
}

BivariateNormal::~BivariateNormal() {

}

double BivariateNormal::pdf(double x, double y) {
        
        return gsl_ran_bivariate_gaussian_pdf(x-mu_x, y-mu_y, sigma_x, sigma_y, rho);
}

double BivariateNormal::pdf(vector<double> x) {
        
        return gsl_ran_bivariate_gaussian_pdf(x[0]-mu_x, x[1]-mu_y, sigma_x, sigma_y, rho);
}

void BivariateNormal::sample(double* x, double* y) {

        gsl_ran_bivariate_gaussian(_r, sigma_x, sigma_y, rho, x, y);
        *x += mu_x;
        *y += mu_y;
}

void BivariateNormal::sample(vector<double>& x) {
        double y1, y2;
        sample(&y1, &y2);
        x[0] = y1; x[1] = y2;
}
