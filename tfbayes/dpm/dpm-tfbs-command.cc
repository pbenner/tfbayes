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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <sstream>

#include <dpm-tfbs-command.hh>
#include <dpm-tfbs-sampler.hh>
#include <dpm-tfbs-state.hh>

using namespace std;

print_cluster_elements_t::print_cluster_elements_t(ssize_t cluster_tag)
        : _cluster_tag(cluster_tag)
{ }

string
print_cluster_elements_t::operator()(const dpm_tfbs_state_t& state, dpm_tfbs_sampler_t& sampler) const {
        stringstream ss;

        for (dpm_tfbs_state_t::const_iterator it = state.begin(); it != state.end(); it++) {
                cluster_t& cluster = **it;
                if (_cluster_tag == cluster.cluster_tag()) {
                        ss << cluster;
                }
        }
        return ss.str();
}
