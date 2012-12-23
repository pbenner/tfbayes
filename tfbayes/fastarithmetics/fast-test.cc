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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <gsl/gsl_sf_gamma.h>

#include <fast-lngamma.hh>
#include <fast-lnbeta.hh>

using namespace std;

int main(void)
{
        // for (size_t i = 1; i <= 100000; i++)
        // {
        //         cout << "       {";
        //         cout.precision(5);
        //         cout.width(8);
        //         cout << 0.01*i
        //              << ", ";
        //         cout.precision(5);
        //         cout.width(8);
        //         cout << gsl_sf_lngamma(0.01*i)
        //              << "},"
        //              << endl;
        // }
        cout << "     lngamma: " << gsl_sf_lngamma(1000) << endl;
        cout << "fast_lngamma: " <<   fast_lngamma(1000) << endl;

        boost::array<double, 5> alpha;
        alpha[0] = 10;
        alpha[1] = 20;
        alpha[2] = 30;
        alpha[3] = 40;
        alpha[4] = 50;

        cout << "fast_lnbeta : " <<   fast_lnbeta<5>(alpha) << endl;
}
