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

        using std::vector<std::vector<double> >::operator=;

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
                // initialize data structures
                likelihood = matrix_t();
                prediction = matrix_t();
                forward    = matrix_t();
                // run forward and backward algorithm
                run_forward (alignment);
                run_backward(alignment);
                // compute posterior marginals
                size_t length = likelihood.size();
                operator=(matrix_t(length-1, vector_t(dim, 0.0)));
                for (size_t k = 1; k < likelihood.size(); k++) {
                        for (size_t i = 0; i < dim; i++) {
                                operator[](k-1)[i] = forward[k][i]*backward[k][i];
                        }
                        normalize(operator[](k-1));
                }
        }

private:
        void run_forward(alignment_t<CODE_TYPE>& alignment) {

                likelihood.push_back(vector_t(dim, 1.0));
                prediction.push_back(vector_t(dim, 0.0));
                // initialize prediction
                for (size_t i = 0; i < dim; i++) {
                        for (size_t j = 0; j < dim; j++) {
                                prediction.back()[i] += transition[j][i]*px_0[j];
                        }
                }
                forward.push_back(px_0);

                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        vector_t likelihood_k(dim, 0.0); // p(z_k | x_k)
                        vector_t prediction_k(dim, 0.0); // p(x_k | z_1:k-1)
                        vector_t forward_k(dim, 0.0);

                        // compute likelihood
                        for (size_t j = 0; j < dim; j++) {
                                likelihood_k[j] = exp(pt_marginal_likelihood(alignment.tree, priors[j]));
                        }

                        // compute prediction
                        for (size_t i = 0; i < dim; i++) {
                                for (size_t j = 0; j < dim; j++) {
                                        prediction_k[i] +=
                                                transition[j][i]*likelihood.back()[j]*prediction.back()[j];
                                }
                        }
                        normalize(prediction_k);

                        // save result
                        likelihood.push_back(likelihood_k);
                        prediction.push_back(prediction_k);

                        // compute p(x_k | z_1:k)
                        for (size_t j = 0; j < dim; j++) {
                                forward_k[j] = likelihood_k[j] * prediction_k[j];
                        }
                        normalize(forward_k);

                        forward.push_back(forward_k);
                }
        }
        void run_backward(alignment_t<CODE_TYPE>& alignment) {
                size_t length = likelihood.size();
                backward = matrix_t(length, vector_t(dim, 0.0));
                backward[length-1] = likelihood[length-1];
                normalize(backward[length-1]);
                // loop through the alignment (backwards)
                for (size_t k = 1; k < length; k++) {
                        vector_t& backward_k   = backward  [length-k-1];
                        vector_t& likelihood_k = likelihood[length-k-1];
                        // loop over states
                        for (size_t i = 0; i < dim; i++) {
                                for (size_t j = 0; j < dim; j++) {
                                        backward_k[i] += likelihood_k[i]*transition[i][j]*backward[length-k][j];
                                }
                        }
                        // normalize to keep numerical precision
                        normalize(backward_k);
                }
        }
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
        matrix_t forward;
        matrix_t backward;
};

#endif /* PHYLOTREE_HMM_HH */
