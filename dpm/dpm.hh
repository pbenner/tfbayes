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

#ifndef DPM_HH
#define DPM_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

using namespace std;

#include "data.hh"
#include "cluster.hh"
#include "statistics.hh"

class DPM {
public:
        // constructors and destructors
        DPM(Data* data);
        virtual ~DPM();

        // operators
        friend ostream& operator<<(std::ostream& o, DPM const& dpm);
              Cluster::cluster& operator[](int c)       { return cl[c]; }
        const Cluster::cluster& operator[](int c) const { return cl[c]; }

        // methods
        bool sample(Data::element& element);
        void gibbsSample(unsigned int steps);
        Cluster::size_type num_clusters() { return cl.size(); }

        virtual Distribution& posteriorPredictive(const Cluster::cluster& cluster) {
                return *posteriorPredictiveDist;
        }
        virtual Distribution& predictive() {
                return *predictiveDist;
        }
        virtual double likelihood() {
                return 0;
        }
        virtual void compute_statistics() {
                hist_likelihood.push_back(likelihood());
                hist_num_clusters.push_back(cl.size());
        }
        virtual Data& get_data() {
                return *da;
        }
        Cluster& get_clusters() {
                return cl;
        }
        vector<double>& get_hist_switches() {
                return hist_switches;
        }
        vector<double>& get_hist_likelihood() {
                return hist_likelihood;
        }

protected:
        Data* da;
        Cluster cl;
        double alpha;

        // distributions
        Distribution* predictiveDist;
        Distribution* posteriorPredictiveDist;

        // gibbs sampler history
        vector<double> hist_switches;
        vector<double> hist_likelihood;
        vector<Cluster::size_type> hist_num_clusters;
};

#endif /* DPM_HH */
