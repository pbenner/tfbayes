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

#ifndef MIXTURE_WEIGHTS_HH
#define MIXTURE_WEIGHTS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <string.h>

class mixture_weights_t : public clonable {
public:
        virtual ~mixture_weights_t() {}

        virtual mixture_weights_t* clone() const = 0;

        virtual const double& operator[](size_t i) const = 0;

        virtual void init(const std::vector<int>& codes) = 0;
        virtual void update(int code, const double *alpha, const double *counts, const double *counts_sum) = 0;
};

class decay_weights_t : public mixture_weights_t {
public:
        virtual decay_weights_t* clone() const {
                return new decay_weights_t(*this);
        }

        virtual const double& operator[](size_t i) const {
                return weights[i];
        }
        virtual void init(const std::vector<int>& codes) {
                const size_t n = codes.size();
                std::vector<double> weights(n, 0.0);
                double k = 1.0;

                for (size_t i = 0; i < n; i++) {
                        if (i < n-1) { k *= 0.5; }
                        weights[i] = k;
                }
                this->weights = weights;
        }
        virtual void update(int code, const double *alpha, const double *counts, const double *counts_sum) {
        }
private:
        std::vector<double> weights;
};

class entropy_weights_t : public mixture_weights_t {
public:
        entropy_weights_t(size_t alphabet_size, size_t max_context, size_t length)
                : alphabet_size(alphabet_size),
                  max_context(max_context),
                  length(length),
                  entropy_max(log(alphabet_size)) {
                entropy = (double*)malloc(length/alphabet_size*sizeof(double));
        }
        entropy_weights_t(const entropy_weights_t& weights)
                : alphabet_size(weights.alphabet_size),
                  max_context(weights.max_context),
                  length(weights.length),
                  entropy_max(weights.entropy_max),
                  weights(weights.weights) {
                entropy = (double*)malloc(length/alphabet_size*sizeof(double));
                memcpy(entropy, weights.entropy, length/alphabet_size*sizeof(double));
                for (size_t i = 0; i < length/alphabet_size; i++) {
                        entropy[i] = 1;
                }
        }
        virtual ~entropy_weights_t() {
                free(entropy);
        }

        virtual entropy_weights_t* clone() const {
                return new entropy_weights_t(*this);
        }

        virtual const double& operator[](size_t i) const {
                return weights[i];
        }
        virtual void init(const std::vector<int>& codes) {
                const size_t n = codes.size();
                double weights_sum = 0.0;
                std::vector<double> weights(n, 0.0);

                for (size_t i = 0; i < n; i++) {
                        weights[i] = (1.0 - entropy[codes[i]/alphabet_size])*(1 - weights_sum);
                        weights_sum += weights[i];
                }
                for (size_t i = 0; i < n; i++) {
                        weights[i] /= weights_sum;
                }
                this->weights = weights;
        }
        virtual void update(int code, const double *alpha, const double *counts, const double *counts_sum) {
                const int from = code - (code%alphabet_size);
                const int k    = code/alphabet_size;

                entropy[k] = 0;
                for (size_t i = from; i < from+alphabet_size; i++) {
                        const double p = (counts[i]+alpha[i])/counts_sum[k];
                        entropy[k] -= p*log(p);
                }
                entropy[k] /= entropy_max;
        }
private:
        size_t  alphabet_size;
        size_t  max_context;
        size_t  length;
        double *entropy;
        double  entropy_max;

        std::vector<double> weights;
};

#endif /* MIXTURE_WEIGHTS_HH */
