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

#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

using namespace std;

extern gsl_rng* _r;

class Distribution {

public:
        Distribution() {}
        virtual ~Distribution() {}

        virtual double pdf(vector<double> x) { return 0.0; }
        virtual double pdf(double x, double y) { return 0.0; }
        virtual void sample(vector<double>& x) {}
        virtual void sample(double* x, double* y) {}
        virtual void update(const gsl_matrix* cov,
                            const gsl_vector* mu) {}
        virtual void update(const vector<double>& mu, 
                            const vector<double>& sigma,
                            double rho) {}
        virtual void update(double mu_x, double mu_y,
                            double sigma_x, double sigma_y,
                            double rho) {}
};

class BivariateNormal : public Distribution {
public:
        BivariateNormal();
        BivariateNormal(const gsl_matrix* cov,
                        const gsl_vector* mu);
        BivariateNormal(const vector<double>& mu, 
                        const vector<double>& sigma,
                        double rho);
        BivariateNormal(double mu_x, double mu_y,
                        double sigma_x, double sigma_y,
                        double rho);
        ~BivariateNormal();

        void update(const gsl_matrix* cov,
                    const gsl_vector* mu);
        void update(const vector<double>& mu, 
                    const vector<double>& sigma,
                    double rho);
        void update(double mu_x, double mu_y,
                    double sigma_x, double sigma_y,
                    double rho);

        double pdf(vector<double> x);
        double pdf(double x, double y);
        void sample(vector<double>& x);
        void sample(double* x, double* y);

private:
        double mu_x;
        double mu_y;
        double sigma_x;
        double sigma_y;
        double rho;
};

#endif /* STATISTICS_HH */
