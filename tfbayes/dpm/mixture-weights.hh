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

        virtual void reset() = 0;
        virtual double next(int code, size_t c, size_t c_max) = 0;
        virtual void update(int code, const double *alpha, const double *counts, const double *counts_sum) = 0;
};

class decay_weights_t : public mixture_weights_t {
public:
        decay_weights_t()
                : weight(1.0) {
        }

        virtual decay_weights_t* clone() const {
                return new decay_weights_t(*this);
        }

        virtual void reset() {
                weight = 1.0;
        }
        virtual double next(int code, size_t c, size_t c_max) {
                if (c != c_max) {
                        weight *= 0.5;
                }
                return weight;
        }
        virtual void update(int code, const double *alpha, const double *counts, const double *counts_sum) {
        }
private:
        double weight;
};

class entropy_weights_t : public mixture_weights_t {
public:
        entropy_weights_t(size_t length, size_t alphabet_size)
                : length(length), alphabet_size(alphabet_size),
                  weight(0.0), weight_sum(0.0), entropy_max(log(alphabet_size)) {
                entropy = (double*)malloc(length/alphabet_size*sizeof(double));
        }
        entropy_weights_t(const entropy_weights_t& weights)
                : length(weights.length),
                  alphabet_size(weights.alphabet_size),
                  weight(weights.weight),
                  weight_sum(weights.weight_sum),
                  entropy_max(weights.entropy_max) {
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

        virtual void reset() {
                weight = 0.0;
                weight_sum = 0.0;
        }
        virtual double next(int code, size_t c, size_t c_max) {
                if (c == c_max) {
                        weight = 1 - weight_sum;
                }
                else {
                        weight = (1 - entropy[code/alphabet_size])*(1 - weight_sum);
                }
                weight_sum += weight;

                return weight;
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
        size_t  length;
        size_t  alphabet_size;
        double  weight;
        double  weight_sum;
        double *entropy;
        double  entropy_max;
};

#endif /* MIXTURE_WEIGHTS_HH */
