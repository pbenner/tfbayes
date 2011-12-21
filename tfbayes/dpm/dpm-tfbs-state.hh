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

#ifndef DPM_TFBS_STATE_HH
#define DPM_TFBS_STATE_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <clonable.hh>

class dpm_tfbs_state_t : public mixture_state_t {
public:
        dpm_tfbs_state_t(const std::vector<size_t> sizes)
                : mixture_state_t(cluster_assignments),
                  cluster_assignments(sizes, -1),
                  // starting positions of tfbs
                  tfbs_start_positions(sizes, 0),
                  // number of transcription factor binding sites
                  num_tfbs(0)
                { }

        // assignments of nucleotides to clusters
        sequence_data_t<cluster_tag_t> cluster_assignments;

        // record start positions of tfbs
        sequence_data_t<short> tfbs_start_positions;

        // posterior distribution
        posterior_t posterior;
        Graph tfbs_graph;

        // keep track of the number of transcription factor binding sites
        size_t num_tfbs;
};

#endif /* DPM_TFBS_STATE_HH */
