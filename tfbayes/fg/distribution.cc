/* Copyright (C) 2013 Philipp Benner
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

#include <tfbayes/fg/distribution.hh>

// normal statistics
////////////////////////////////////////////////////////////////////////////////

exponential_family_i::vector_t
normal_distribution_t::statistics(double x)  {
        vector_t T(2, 0.0);
        T[0] += x;
        T[1] += x*x;
        return T;
}

exponential_family_i::vector_t
normal_distribution_t::statistics(const vector_t& x) {
        return statistics(x[0]);
}

// gamma statistics
////////////////////////////////////////////////////////////////////////////////

exponential_family_i::vector_t
gamma_distribution_t::statistics(double x) {
        vector_t T(2, 0.0);
        T[0] += std::log(x);
        T[1] += x;
        return T;
}

exponential_family_i::vector_t
gamma_distribution_t::statistics(const vector_t& x) {
        return statistics(x[0]);
}

// dirichlet
////////////////////////////////////////////////////////////////////////////////

exponential_family_i::vector_t
dirichlet_distribution_t::statistics(const vector_t& x) {
        vector_t T(x.size(), 0.0);
        for (size_t i = 0; i < x.size(); i++) {
                T[i] = std::log(x[i]);
        }
        return T;
}

// discrete
////////////////////////////////////////////////////////////////////////////////

exponential_family_i::vector_t
categorical_distribution_t::statistics(const vector_t& x) {
        return x;
}
