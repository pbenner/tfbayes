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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <math.h>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

#include <data.hh>
#include <distribution.hh>

using namespace std;

gsl_rng* _r;

ProductDirichlet::ProductDirichlet(gsl_matrix* _alpha)
{
        for (size_t i = 0; i < _alpha->size1; i++) {
                double sum = 0;
                alpha.push_back (vector<double>(_alpha->size2+1, 0));
                counts.push_back(vector<size_t>(_alpha->size2+1, 0));
                for (size_t j = 0; j < _alpha->size2; j++) {
                        alpha[i][j] = gsl_matrix_get(_alpha, i, j);
                        sum += alpha[i][j];
                }
                alpha[i][_alpha->size2] = sum;
        }
}

ProductDirichlet::ProductDirichlet(const ProductDirichlet& distribution)
        : alpha(distribution.alpha), counts(distribution.counts)
{
}

ProductDirichlet::~ProductDirichlet() {
}

ProductDirichlet*
ProductDirichlet::clone() const {
        return new ProductDirichlet(*this);
}

size_t
ProductDirichlet::count_observations(const word_t& word) const {
        return word.length/counts.size();
}

size_t
ProductDirichlet::remove_observations(const word_t& word) {
        for (size_t i = 0; i < word.length; i += counts.size()) {
                for (size_t j = 0; j < counts.size(); j++) {
                        switch (word.sequences[word.sequence][word.position+i+j]) {
                        case 'A':
                        case 'a':
                                counts[j][0]--;
                        break;
                        case 'C':
                        case 'c':
                                counts[j][1]--;
                        break;
                        case 'G':
                        case 'g':
                                counts[j][2]--;
                        break;
                        case 'T':
                        case 't':
                                counts[j][3]--;
                        break;
                        }
                        counts[j][4]--;
                }
        }
        return word.length/counts.size();
}

size_t
ProductDirichlet::add_observations(const word_t& word) {
        for (size_t i = 0; i < word.length; i += counts.size()) {
                for (size_t j = 0; j < counts.size(); j++) {
                        switch (word.sequences[word.sequence][word.position+i+j]) {
                        case 'A':
                        case 'a':
                                counts[j][0]++;
                        break;
                        case 'C':
                        case 'c':
                                counts[j][1]++;
                        break;
                        case 'G':
                        case 'g':
                                counts[j][2]++;
                        break;
                        case 'T':
                        case 't':
                                counts[j][3]++;
                        break;
                        }
                        counts[j][4]++;
                }
        }
        return word.length/counts.size();
}

double ProductDirichlet::pdf(const word_t& word) const {
        double result = 1;

        for (size_t i = 0; i < word.length; i += counts.size()) {
                for (size_t j = 0; j < counts.size(); j++) {
                        switch (word.sequences[word.sequence][word.position+i+j]) {
                        case 'A':
                        case 'a':
                                result *= (counts[j][0]+alpha[j][0])/(counts[j][4]+alpha[j][4]);
                        break;
                        case 'C':
                        case 'c':
                                result *= (counts[j][1]+alpha[j][1])/(counts[j][4]+alpha[j][4]);
                        break;
                        case 'G':
                        case 'g':
                                result *= (counts[j][2]+alpha[j][2])/(counts[j][4]+alpha[j][4]);
                        break;
                        case 'T':
                        case 't':
                                result *= (counts[j][3]+alpha[j][3])/(counts[j][4]+alpha[j][4]);
                        break;
                        }
                }
        }
        return result;
}

double ProductDirichlet::log_pdf(const word_t& word) const {
        double result = 0;

        for (size_t i = 0; i < word.length; i += counts.size()) {
                for (size_t j = 0; j < counts.size(); j++) {
                        switch (word.sequences[word.sequence][word.position+i+j]) {
                        case 'A':
                        case 'a':
                                result += logl((counts[j][0]+alpha[j][0])/(counts[j][4]+alpha[j][4]));
                        break;
                        case 'C':
                        case 'c':
                                result += logl((counts[j][1]+alpha[j][1])/(counts[j][4]+alpha[j][4]));
                        break;
                        case 'G':
                        case 'g':
                                result += logl((counts[j][2]+alpha[j][2])/(counts[j][4]+alpha[j][4]));
                        break;
                        case 'T':
                        case 't':
                                result += logl((counts[j][3]+alpha[j][3])/(counts[j][4]+alpha[j][4]));
                        break;
                        }
                }
        }
        return result;
}

//
// \sum_{x \in X} n_x log(\frac{n_x + \alpha_x}{\sum_{x' \in X} n_{x'} + \alpha_{x'}})
//
double ProductDirichlet::log_likelihood() const {
        double result = 0;

        for (size_t j = 0; j < counts.size(); j++) {
                for (size_t k = 0; k < counts[j].size(); k++) {
                        result += counts[j][k]*logl((counts[j][k]+alpha[j][k])/(counts[j][4]+alpha[j][4]));
                }
        }
        return result;
}
