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

#ifndef __TFBAYES_DPM_MIXTURE_WEIGHTS_HH__
#define __TFBAYES_DPM_MIXTURE_WEIGHTS_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <cmath>

#include <tfbayes/utility/clonable.hh>

class mixture_weights_t : virtual public std::vector<double>, public virtual clonable {
public:
        virtual ~mixture_weights_t() {}

        virtual mixture_weights_t* clone() const = 0;

        virtual void init(const std::vector<int>& codes) = 0;
        virtual void update(int code, const double *alpha, const double *counts, const double *counts_sum) = 0;
};

class decay_weights_t : public mixture_weights_t {
public:
        decay_weights_t(size_t max_context)
                : std::vector<double>(max_context+1, 0.0) { }

        virtual decay_weights_t* clone() const {
                return new decay_weights_t(*this);
        }

        virtual void init(const std::vector<int>& codes) {
                const size_t n = codes.size();
                double k = 1.0;

                for (size_t i = 0; i < n; i++) {
                        if (i < n-1) { k *= 0.5; }
                        operator[](i) = k;
                }
        }

        virtual void update(int code, const double *alpha, const double *counts, const double *counts_sum) {
        }
};

class entropy_weights_t : public mixture_weights_t {
public:
        entropy_weights_t(size_t alphabet_size, size_t max_context, size_t length)
                : std::vector<double>(max_context+1, 0.0),
                  alphabet_size(alphabet_size),
                  max_context(max_context),
                  length(length),
                  entropy_max(log(alphabet_size)) {
                entropy = (double*)malloc(length/alphabet_size*sizeof(double));
                for (size_t i = 0; i < length/alphabet_size; i++) {
                        entropy[i] = 1.0;
                }
        }
        entropy_weights_t(const entropy_weights_t& weights)
                : std::vector<double>(weights),
                  alphabet_size(weights.alphabet_size),
                  max_context(weights.max_context),
                  length(weights.length),
                  entropy_max(weights.entropy_max) {
                entropy = (double*)malloc(length/alphabet_size*sizeof(double));
                for (size_t i = 0; i < length/alphabet_size; i++) {
                        entropy[i] = weights.entropy[i];
                }
        }
        virtual ~entropy_weights_t() {
                free(entropy);
        }

        virtual entropy_weights_t& operator=(const entropy_weights_t& weights) {
                std::vector<double>::operator=(weights);
                for (size_t i = 0; i < length/alphabet_size; i++) {
                        entropy[i] = weights.entropy[i];
                }
                return *this;
        }

        virtual entropy_weights_t* clone() const {
                return new entropy_weights_t(*this);
        }

        virtual void init(const std::vector<int>& codes) {
                const size_t n = codes.size();
                double certainty_sum1 = 0.0;
                double certainty_sum2 = 0.0;
                double weights_sum = 0.0;

                for (size_t i = 0; i < n; i++) {
                        const double h  = entropy[codes[i]/alphabet_size];
                        certainty_sum1 += (2.0 - h);
                }
                for (size_t i = 0; i < n; i++) {
                        const double h  = entropy[codes[i]/alphabet_size];
                        const double c  = (2.0 - h)/certainty_sum1;
                        operator[](i)   = c*(1.0 - certainty_sum2);
                        weights_sum    += operator[](i);
                        certainty_sum2 += c;
                }
                for (size_t i = 0; i < n; i++) {
                        operator[](i) /= weights_sum;
                }
        }

        virtual void update(int code, const double *alpha, const double *counts, const double *counts_sum) {
                const int from = code - (code%alphabet_size);
                const int k    = code/alphabet_size;

                entropy[k] = 0.0;
                for (size_t i = from; i < from+alphabet_size; i++) {
                        const double p = (counts[i]+alpha[i])/counts_sum[k];
                        entropy[k] -= p*log(p);
                }
                entropy[k] /= entropy_max;
        }
private:
        const size_t  alphabet_size;
        const size_t  max_context;
        const size_t  length;
        const double  entropy_max;
        double *entropy;
};

#endif /* __TFBAYES_DPM_MIXTURE_WEIGHTS_HH__ */
