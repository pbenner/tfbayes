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

#ifndef __TFBAYES_ENTROPY_ENTROPY_HH__
#define __TFBAYES_ENTROPY_ENTROPY_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cmath>

#include <boost/foreach.hpp>

template <class input_type = double>
double entropy(const std::vector<input_type>& v) {
        double result = 0.0;
        BOOST_FOREACH(const input_type& i, v) {
                result -= static_cast<double>(i)*std::log(i);
        }
        return result;
}


#endif /* __TFBAYES_ENTROPY_ENTROPY_HH__ */
