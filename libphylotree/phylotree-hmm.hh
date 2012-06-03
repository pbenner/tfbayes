/* Copyright (C) 2012 Philipp Benner
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

#ifndef PHYLOTREE_HMM_HH
#define PHYLOTREE_HMM_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <math.h>

#include <alignment.hh>
#include <phylotree.hh>
#include <phylotree-parser.hh>
#include <phylotree-polynomial.hh>
#include <marginal-likelihood.hh>

#include <tfbayes/exception.h>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class phylotree_hmm_t : public std::vector<std::vector<double> >
{
public:
        typedef std::vector<double>   vector_t;
        typedef std::vector<vector_t> matrix_t;
        typedef std::vector<exponent_t<CODE_TYPE, ALPHABET_SIZE> > priors_t;

        phylotree_hmm_t(
                const vector_t& px_0,
                const matrix_t& transition,
                const priors_t& priors)
                : std::vector<std::vector<double> >(),
                  px_0(px_0),
                  transition(transition),
                  priors(priors) {

                dim = px_0.size();
        }

        void run(alignment_t<CODE_TYPE>& alignment) {

                likelihood.push_back(vector_t(dim, 1.0));
                prediction.push_back(vector_t(dim, 0.0));
                // initialize prediction
                for (size_t i = 0; i < dim; i++) {
                        for (size_t j = 0; j < dim; j++) {
                                prediction.back()[i] += transition[j][i]*px_0[j];
                        }
                }

                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        vector_t likelihood_i(dim, 0.0); // p(z_k | x_k)
                        vector_t prediction_i(dim, 0.0); // p(x_k | z_1:k-1)
                        vector_t p_hidden(dim, 0.0);

                        // compute likelihood
                        for (size_t j = 0; j < dim; j++) {
                                likelihood_i[j] = exp(pt_marginal_likelihood(alignment.tree, priors[j]));
                        }

                        // compute prediction
                        for (size_t i = 0; i < dim; i++) {
                                for (size_t j = 0; j < dim; j++) {
                                        prediction_i[i] +=
                                                transition[j][i]*likelihood.back()[j]*prediction.back()[j];
                                }
                        }
                        normalize(prediction_i);

                        // save result
                        likelihood.push_back(likelihood_i);
                        prediction.push_back(prediction_i);

                        // compute p(x_k | z_1:k)
                        for (size_t j = 0; j < dim; j++) {
                                p_hidden[j] = likelihood_i[j] * prediction_i[j];
                        }
                        normalize(p_hidden);

                        push_back(p_hidden);
                }
        }

private:
        void normalize(vector_t& vector) const {
                double normalization = 0.0;

                for (size_t i = 0; i < vector.size(); i++) {
                        normalization += vector[i];
                }
                for (size_t i = 0; i < vector.size(); i++) {
                        vector[i] /= normalization;
                }
        }

        size_t dim;

        vector_t px_0;
        matrix_t transition;
        priors_t priors;

        matrix_t likelihood;
        matrix_t prediction;
};

#endif /* PHYLOTREE_HMM_HH */
