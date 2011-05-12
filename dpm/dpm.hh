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

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include "data.hh"
#include "dpm.hh"
#include "clusters.hh"
#include "statistics.hh"

using namespace std;

class DPM {
public:
        DPM(Data* data);
        ~DPM();

        double likelihood();
        void compute_statistics();

        Data& get_data() {
                return *static_cast<Data *>(da);
        }

        bool sample(Data::element& element);
        void gibbsSample(unsigned int steps);

        // operators
        ////////////////////////////////////////////////////////////////////////
        friend ostream& operator<<(std::ostream& o, DPM const& dpm);
              Clusters::cluster& operator[](int c)       { return cl[c]; }
        const Clusters::cluster& operator[](int c) const { return cl[c]; }

        // methods
        ////////////////////////////////////////////////////////////////////////
        Clusters::size_type num_clusters() {
                return cl.size();
        }
        Clusters& get_clusters() {
                return cl;
        }
        vector<double>& get_hist_switches() {
                return hist_switches;
        }
        vector<double>& get_hist_likelihood() {
                return hist_likelihood;
        }
        gsl_matrix* get_posterior() {
                return posterior;
        }

        // constants
        ////////////////////////////////////////////////////////////////////////
        static const int BG_CLUSTER  = 0;
        static const int BG_LENGTH   = 1;
        static const int TFBS_LENGTH = 10;
        static const int NUCLEOTIDES = 4;

private:
        // data and clusters
        Data* da;
        Clusters cl;

        // parameters
        double alpha;
        double lambda;

        // priors
        gsl_matrix* tfbs_alpha;
        gsl_matrix* bg_alpha;

        // gibbs sampler history
        unsigned int total_sampling_steps;
        vector<double> hist_switches;
        vector<double> hist_likelihood;
        vector<Clusters::size_type> hist_num_clusters;
        gsl_matrix* posterior;

        // private methods
        void count_statistic(const Clusters::cluster& cluster, gsl_matrix* alpha, gsl_matrix* counts);
        bool check_element(Data::element& element);
        void assign_block(char* nucleotides, Data::element& element, Clusters::cluster& c);
        void release_block(char* nucleotides, Data::element& element, Clusters::cluster& c);
        int num_tfbs();
};

#endif /* DPM_HH */
