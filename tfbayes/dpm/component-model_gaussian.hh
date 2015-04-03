/* Copyright (C) 2011-2013 Philipp Benner
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

#ifndef __TFBAYES_DPM_COMPONENT_MODEL_GAUSSIAN_HH__
#define __TFBAYES_DPM_COMPONENT_MODEL_GAUSSIAN_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

// Bivariate Gaussian
////////////////////////////////////////////////////////////////////////////////

class bivariate_normal_t : public component_model_t {
public:
         bivariate_normal_t();
         bivariate_normal_t(const std::matrix<double>& Sigma,
                            const std::matrix<double>& Sigma_0,
                            const std::vector<double>& mu_0,
                            const data_t<std::vector<double> >& data);
         bivariate_normal_t(const bivariate_normal_t& bn);
        ~bivariate_normal_t();

        bivariate_normal_t* clone() const;

        friend void swap(bivariate_normal_t& first, bivariate_normal_t& second) {
                using std::swap;
                swap(static_cast<component_model_t&>(first),
                     static_cast<component_model_t&>(second));
                swap(first._Sigma_0_inv, second._Sigma_0_inv);
                swap(first._mu_0,        second._mu_0);
                swap(first._Sigma,       second._Sigma);
                swap(first._Sigma_inv,   second._Sigma_inv);
                swap(first._mu,          second._mu);
                swap(first._N,           second._N);
                swap(first._Sigma_N,     second._Sigma_N);
                swap(first._Sigma_N_inv, second._Sigma_N_inv);
                swap(first._mu_N,        second._mu_N);
                swap(first._dimension,   second._dimension);
                swap(first._inv_perm,    second._inv_perm);
                swap(first._inv_tmp,     second._inv_tmp);
                swap(first._tmp1,        second._tmp1);
                swap(first._tmp2,        second._tmp2);
                swap(first._data,        second._data);
        }

        bivariate_normal_t& operator=(const component_model_t& component_model);

        size_t add(const range_t& range);
        size_t remove(const range_t& range);
        size_t count(const range_t& range);
        double predictive(const range_t& range);
        double predictive(const std::vector<range_t>& range_set);
        double log_predictive(const range_t& range);
        double log_predictive(const std::vector<range_t>& range_set);
        double log_likelihood() const;
        std::vector<double> mean() const;

        const data_t<std::vector<double> >& data() const {
                return *_data;
        }

        friend std::ostream& operator<< (std::ostream& o, const bivariate_normal_t& pd);

protected:
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
        size_t _dimension;
        gsl_permutation* _inv_perm;
        gsl_matrix* _inv_tmp;
        gsl_vector* _tmp1;
        gsl_vector* _tmp2;

        void inverse(gsl_matrix* dst, const gsl_matrix* src);
        void update();

        const data_t<std::vector<double> >* _data;
};

#endif /* __TFBAYES_DPM_COMPONENT_MODEL_GAUSSIAN_HH__ */
