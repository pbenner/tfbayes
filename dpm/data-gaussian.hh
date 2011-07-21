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

#ifndef DATA_GAUSSIAN_HH
#define DATA_GAUSSIAN_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <datatypes.hh>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Data_Gaussian : public data_t<std::vector<double> > {
public:
         Data_Gaussian(size_t cluster, size_t samples, gsl_matrix* Sigma);
        ~Data_Gaussian();

        size_t elements() const;
        size_t length() const;

private:
        const size_t _elements;
        const size_t _length;
        const size_t _cluster;

        gsl_vector** _mu;

        std::vector<std::vector<double> > generate_samples(size_t cluster, size_t samples, gsl_matrix* Sigma);
};

#endif /* DATA_GAUSSIAN_HH */
