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

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#include <clonable.hh>
#include <datatypes.hh>

extern gsl_rng* _r;

class Distribution : public clonable {

public:
        Distribution() {}
        Distribution(const Distribution& distribution) {
                std::cout << "Distribution copy constructor called." << std::endl;
                exit(EXIT_FAILURE);
        }

        // datatypes
        typedef std::vector<std::vector<double> > alpha_t;
        typedef std::vector<std::vector<size_t> > counts_t;

        // purely virtual functions
        virtual size_t add_observations(const word_t& word) = 0;
        virtual size_t remove_observations(const word_t& word) = 0;
        virtual size_t count_observations(const word_t& word) const = 0;
        virtual double pdf(const word_t& word) const = 0;
        virtual double log_pdf(const word_t& word) const = 0;
        virtual double log_likelihood() const = 0;

        virtual Distribution* clone() const = 0;
};

class ProductDirichlet : public Distribution {
public:
         ProductDirichlet(gsl_matrix* alpha);
         ProductDirichlet(const ProductDirichlet& distribution);
        ~ProductDirichlet();

        size_t add_observations(const word_t& word);
        size_t remove_observations(const word_t& word);
        size_t count_observations(const word_t& word) const;
        double pdf(const word_t& word) const;
        double log_pdf(const word_t& word) const;
        double log_likelihood() const;

        ProductDirichlet* clone() const;

private:
        alpha_t alpha;
        counts_t counts;
};

#endif /* STATISTICS_HH */
