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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cmath>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/linalg.hh>

template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class phylotree_hmm_t : public std::matrix<double>
{
public:
        typedef std::vector<exponent_t<AS, PC> > priors_t;

        using std::matrix<double>::operator=;

        phylotree_hmm_t(
                const std::vector<double>& px_0,
                const std::matrix<double>& transition,
                const priors_t& priors)
                : std::matrix<double>(),
                  px_0(px_0),
                  transition(transition),
                  priors(priors) {

                dim = px_0.size();
        }

        void run(const pt_root_t& tree, const alignment_t<AC>& alignment) {
                // initialize data structures
                likelihood = std::matrix<double>();
                prediction = std::matrix<double>();
                forward    = std::matrix<double>();
                // run forward and backward algorithm
                run_forward (tree, alignment);
                run_backward();
                // compute posterior marginals
                size_t length = likelihood.size();
                operator=(std::matrix<double>(length-1, dim));
                for (size_t k = 1; k < likelihood.size(); k++) {
                        for (size_t i = 0; i < dim; i++) {
                                operator[](k-1)[i] = forward[k][i]*backward[k][i];
                        }
                        normalize(operator[](k-1));
                }
        }

protected:
        void run_forward(const pt_root_t& tree, const alignment_t<AC>& alignment) {

                likelihood.push_back(std::vector<double>(dim, 1.0));
                prediction.push_back(std::vector<double>(dim, 0.0));
                // initialize prediction
                for (size_t i = 0; i < dim; i++) {
                        for (size_t j = 0; j < dim; j++) {
                                prediction.back()[i] += transition[j][i]*px_0[j];
                        }
                }
                forward.push_back(px_0);

                for (typename alignment_t<AC>::const_iterator it = alignment.begin(); it != alignment.end(); it++) {
                        std::vector<double> likelihood_k(dim, 0.0); // p(z_k | x_k)
                        std::vector<double> prediction_k(dim, 0.0); // p(x_k | z_1:k-1)
                        std::vector<double> forward_k(dim, 0.0);

                        // compute likelihood
                        for (size_t j = 0; j < dim; j++) {
                                likelihood_k[j] = exp(pt_marginal_likelihood(tree, *it, priors[j]));
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
        void run_backward() {
                size_t length = likelihood.size();
                backward = std::matrix<double>(length, dim);
                backward[length-1] = likelihood[length-1];
                normalize(backward[length-1]);
                // loop backwards (we don't need the alignment here)
                for (size_t k = 1; k < length; k++) {
                        std::vector<double>& backward_k   = backward  [length-k-1];
                        std::vector<double>& likelihood_k = likelihood[length-k-1];
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
        void normalize(std::vector<double>& vector) const {
                double normalization = 0.0;

                for (size_t i = 0; i < vector.size(); i++) {
                        normalization += vector[i];
                }
                for (size_t i = 0; i < vector.size(); i++) {
                        vector[i] /= normalization;
                }
        }

        size_t dim;

        std::vector<double> px_0;
        std::matrix<double> transition;
        priors_t priors;

        std::matrix<double> likelihood;
        std::matrix<double> prediction;
        std::matrix<double> forward;
        std::matrix<double> backward;
};

#endif /* PHYLOTREE_HMM_HH */
