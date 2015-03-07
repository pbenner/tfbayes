/* Copyright (C) 2011-2015 Philipp Benner
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

#include <tfbayes/dpm/dpm-tfbs-options.hh>

using namespace std;

ostream&
operator<<(ostream& o, const tfbs_options_t& options) {
        o << "Options:"                   << endl
          << "-> alpha                = " << options.alpha                << endl
          << "-> discount             = " << options.discount             << endl
          << "-> lambda               = " << options.lambda               << endl
          << "-> block samples        = " << options.block_samples        << endl
          << "-> block samples period = " << options.block_samples_period << endl
          << "-> optimize             = " << options.optimize             << endl
          << "-> optimize period      = " << options.optimize_period      << endl
          << "-> initial temperature  = " << options.initial_temperature  << endl
          << "-> process prior        = " << options.process_prior        << endl
          << "-> background model     = " << options.background_model     << endl
          << "-> background context   = " << options.background_context   << endl
          << "-> population_size      = " << options.population_size      << endl
          << "-> socket_file          = " << options.socket_file          << endl
          << "-> verbose              = " << options.verbose              << endl;
        return o;
}
