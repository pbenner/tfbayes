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

#ifndef _PERMUTATION_H_
#define _PERMUTATION_H_

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cstdlib>
#include <cstddef>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

class random_permutation_t : std::vector<double> {
public:
        random_permutation_t(size_t length, boost::random::mt19937& gen)
                : std::vector<double>(length-1, 0.0) {

                boost::uniform_01<> uniform;

                for (size_t i = 0; i < length-1; i++) {
                        operator[](i) = uniform(gen);
                }
        }

        ptrdiff_t operator() (ptrdiff_t max) const {
                return static_cast<ptrdiff_t>(operator[](max-2) * max);
        }
};

#endif /* _PERMUTATION_H_ */
