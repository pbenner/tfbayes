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
#include <data.hh>
#include <datatypes.hh>
#include <nucleotide-sequence.hh>

#include <parsmm/abstract_set.h>
#include <parsmm/static_pars_tree.h>

class ComponentModel : public clonable {
public:
        ComponentModel() {}
        ComponentModel(const ComponentModel& model) {
                std::cout << "Distribution copy constructor called." << std::endl;
                exit(EXIT_FAILURE);
        }

        // purely virtual functions
        virtual size_t add(const range_t& range) = 0;
        virtual size_t remove(const range_t& range) = 0;
        virtual size_t count(const range_t& range) = 0;
        virtual double predictive(const range_t& range) = 0;
        virtual double log_predictive(const range_t& range) = 0;
        virtual double log_likelihood() const = 0;

        virtual ComponentModel* clone() const = 0;
};

class ProductDirichlet : public ComponentModel {
public:
         ProductDirichlet(const std::matrix<double>& alpha, const sequence_data_t<short>& data);
         ProductDirichlet(const ProductDirichlet& distribution);
        ~ProductDirichlet();

        // datatypes
        typedef std::vector<std::vector<size_t> > alpha_t;
        typedef std::vector<std::vector<size_t> > counts_t;

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double log_predictive(const range_t& range);
        double log_likelihood() const;

        ProductDirichlet* clone() const;

        friend std::ostream& operator<< (std::ostream& o, const ProductDirichlet& pd);

private:
        alpha_t alpha;
        counts_t counts;

        const sequence_data_t<short>& _data;

        const size_t _size1;
        const size_t _size2;
};

class ParsimoniousTree : public ComponentModel {
public:
         ParsimoniousTree(size_t alphabet_size, size_t tree_depth,
                          const sequence_data_t<short>& data,
                          const sequence_data_t<cluster_tag_t>& cluster_assignments,
                          cluster_tag_t cluster_tag);
         ParsimoniousTree(const ParsimoniousTree& distribution);
        ~ParsimoniousTree();

        // datatypes
        typedef std::vector<std::vector<size_t> > alpha_t;
        typedef std::vector<std::vector<size_t> > counts_t;

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double log_predictive(const range_t& range);
        double log_likelihood() const;
        ParsimoniousTree* clone() const;

        friend std::ostream& operator<< (std::ostream& o, const ParsimoniousTree& pd);

private:
        abstract_set_t* _as;
        static_pars_tree_t* _pt;
        size_t _counts_length;
        count_t* _counts;

        const sequence_data_t<short>&         _data;
              sequence_data_t<context_t>      _context;
        const sequence_data_t<cluster_tag_t>& _cluster_assignments;

        const cluster_tag_t _cluster_tag;
        const size_t        _tree_depth;

        // internal methods
        size_t max_from_context(const range_t& range) const;
        size_t max_to_context(const range_t& range) const;
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
        double predictive(const range_t& range);
        double log_predictive(const range_t& range);
        double log_likelihood() const;
        const gsl_vector* mean() const;

        BivariateNormal* clone() const;

        friend std::ostream& operator<< (std::ostream& o, const BivariateNormal& pd);

private:
        // prior
        gsl_matrix* _Sigma_0_inv;
        gsl_vector* _mu_0;

        // likelihood
        gsl_matrix* _Sigma;
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
