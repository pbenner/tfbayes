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

#ifndef __TFBAYES_PHYLOTREE_UTILITY_HH__
#define __TFBAYES_PHYLOTREE_UTILITY_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/utility/polynomial.hh>

template<typename T>
std::ostream& operator<< (std::ostream& o, const exponent_t<4, T>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}

template<typename T>
std::ostream& operator<< (std::ostream& o, const exponent_t<5, T>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];
        if(exponent[4]) o << " Pn^" << exponent[4];

        return o;
}

#endif /* __TFBAYES_PHYLOTREE_UTILITY_HH__ */
