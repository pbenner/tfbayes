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

#ifndef COMPONENT_MODEL_HH
#define COMPONENT_MODEL_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>

#include <clonable.hh>
#include <datatypes.hh>

extern gsl_rng* _r;

class ComponentModel : public clonable {
public:
        ComponentModel() {}
        ComponentModel(const ComponentModel& model) {
                std::cout << "Distribution copy constructor called." << std::endl;
                exit(EXIT_FAILURE);
        }

        // datatypes
        typedef std::vector<std::vector<size_t> > alpha_t;
        typedef std::vector<std::vector<size_t> > counts_t;

        // purely virtual functions
        virtual size_t add(const range_t& range) = 0;
        virtual size_t remove(const range_t& range) = 0;
        virtual size_t count(const range_t& range) = 0;
        virtual double pdf(const range_t& range) const = 0;
        virtual double log_pdf(const range_t& range) const = 0;
        virtual double log_likelihood() const = 0;

        virtual ComponentModel* clone() const = 0;
};

class ProductDirichlet : public ComponentModel {
public:
         ProductDirichlet(gsl_matrix* alpha, const sequence_data_t<short>& data);
         ProductDirichlet(const ProductDirichlet& distribution);
        ~ProductDirichlet();

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double pdf(const range_t& range) const;
        double log_pdf(const range_t& range) const;
        double log_likelihood() const;

        ProductDirichlet* clone() const;

private:
        alpha_t alpha;
        counts_t counts;

        const sequence_data_t<short>& _data;

        const size_t _size1;
        const size_t _size2;
};

class BivariateNormal : public ComponentModel {
public:
         BivariateNormal();
         BivariateNormal(const gsl_matrix* Sigma,
                         const gsl_matrix* Sigma_0,
                         const gsl_vector* mu_0,
                         const data_t<std::vector<double> >& data);
         BivariateNormal(const BivariateNormal& bn);
        ~BivariateNormal();

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double pdf(const range_t& range) const;
        double log_pdf(const range_t& range) const;
        double log_likelihood() const;
        const gsl_vector* mean() const;

        BivariateNormal* clone() const;

private:
        // prior
        gsl_matrix* _Sigma_0_inv;
        gsl_vector* _mu_0;

        // likelihood
        gsl_matrix* _Sigma_inv;
        gsl_vector* _mu;
        double _N;

        // posterior
        gsl_matrix* _Sigma_N;
        gsl_matrix* _Sigma_N_inv;
        gsl_vector* _mu_N;

        // other
        const size_t _dimension;
        gsl_permutation* _inv_perm;
        gsl_matrix* _inv_tmp;
        gsl_vector* _tmp1;
        gsl_vector* _tmp2;

        void inverse(gsl_matrix* dst, const gsl_matrix* src);
        void update();

        const data_t<std::vector<double> >& _data;
};

#endif /* COMPONENT_MODEL_HH */
