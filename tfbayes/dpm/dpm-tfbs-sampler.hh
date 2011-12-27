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

#ifndef DPM_TFBS_SAMPLER_HH
#define DPM_TFBS_SAMPLER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <sstream>
#include <vector>

#include <dpm-tfbs.hh>
#include <pmcmc.hh>

class dpm_tfbs_sampler_t : public PopulationMCMC {
public:
        dpm_tfbs_sampler_t(const tfbs_options_t& options, const std::vector<std::string>& sequences, size_t n)
                : PopulationMCMC(n),
                  _data(sequences, options.tfbs_length),
                  _gdpm(n, NULL) {
                for (size_t i = 0; i < _size; i++) {
                        std::stringstream ss; ss << "Sampler " << i+1;
                        _gdpm[i]       = new DpmTfbs(options, _data);
                        _population[i] = new HybridSampler(*_gdpm[i], _gdpm[i]->state(), _data, ss.str(),
                                                           options.metropolis_optimize);
                }
        }

        data_tfbs_t _data;
        std::vector<DpmTfbs*> _gdpm;
};

#endif /* DPM_TFBS_HH */
