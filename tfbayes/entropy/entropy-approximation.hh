/* Copyright (C) 2015 Philipp Benner
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

#ifndef __TFBAYES_ENTROPY_APPROXIMATION_HH__
#define __TFBAYES_ENTROPY_APPROXIMATION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <array>
#include <vector>
#include <cmath>

#include <tfbayes/utility/histogram.hh>

#include <tfbayes/entropy/entropy-approximation-2.hh>
#include <tfbayes/entropy/entropy-approximation-3.hh>
#include <tfbayes/entropy/entropy-approximation-4.hh>
#include <tfbayes/entropy/entropy-approximation-5.hh>
#include <tfbayes/entropy/entropy-approximation-6.hh>
#include <tfbayes/entropy/entropy-approximation-7.hh>
#include <tfbayes/entropy/entropy-approximation-8.hh>
#include <tfbayes/entropy/entropy-approximation-9.hh>
#include <tfbayes/entropy/entropy-approximation-10.hh>

template <class input_type, class result_type>
histogram_t<input_type, result_type> entropy_approximation(size_t k)
{
        switch(k) {
        case  2: return histogram_t<input_type, result_type>(0.0, std::log(k), entropy_histogram_2,  entropy_histogram_2_counts);
        case  3: return histogram_t<input_type, result_type>(0.0, std::log(k), entropy_histogram_3,  entropy_histogram_3_counts);
        case  4: return histogram_t<input_type, result_type>(0.0, std::log(k), entropy_histogram_4,  entropy_histogram_4_counts);
        case  5: return histogram_t<input_type, result_type>(0.0, std::log(k), entropy_histogram_5,  entropy_histogram_5_counts);
        case  6: return histogram_t<input_type, result_type>(0.0, std::log(k), entropy_histogram_6,  entropy_histogram_6_counts);
        case  7: return histogram_t<input_type, result_type>(0.0, std::log(k), entropy_histogram_7,  entropy_histogram_7_counts);
        case  8: return histogram_t<input_type, result_type>(0.0, std::log(k), entropy_histogram_8,  entropy_histogram_8_counts);
        case  9: return histogram_t<input_type, result_type>(0.0, std::log(k), entropy_histogram_9,  entropy_histogram_9_counts);
        case 10: return histogram_t<input_type, result_type>(0.0, std::log(k), entropy_histogram_10, entropy_histogram_10_counts);
        default: throw std::runtime_error("Invalid cardinality!");
        }
}

#endif /* __TFBAYES_ENTROPY_APPROXIMATION_HH__ */
