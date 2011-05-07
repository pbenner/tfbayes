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

#include <data.hh>

using namespace std;

extern gsl_rng* _r;

class Distribution {

public:
        Distribution() {}
        virtual ~Distribution() {}

        virtual double ln_pdf(Data::x_t x) { return 0.0; }
};

class ProductDirichlet : public Distribution {
public:
        ProductDirichlet();
        ProductDirichlet(double lambda, gsl_matrix* counts);
        ~ProductDirichlet();

        void update(double lambda, gsl_matrix* counts);
        double ln_pdf(Data::x_t x);

private:
        double lambda;
        gsl_matrix* alpha;
};

#endif /* STATISTICS_HH */
