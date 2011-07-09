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

#include <clustermanager.hh>
#include <distribution.hh>

class DPM {
public:
         DPM(size_t n, char *sequences[]);
        ~DPM();

        // operators
        ////////////////////////////////////////////////////////////////////////
              Cluster& operator[](cluster_tag_t c)       { return (*cluster_manager)[c]; }
        const Cluster& operator[](cluster_tag_t c) const { return (*cluster_manager)[c]; }

        friend std::ostream& operator<<(std::ostream& o, const DPM& dpm);

        // methods
        ////////////////////////////////////////////////////////////////////////
        std::vector<double>& get_hist_switches() {
                return hist_switches;
        }
        std::vector<double>& get_hist_likelihood() {
                return hist_likelihood;
        }
        const std::vector<std::vector<double> >& get_posterior() const {
                return posterior;
        }

        double compute_likelihood();
        void update_posterior();

        bool valid_for_sampling(const element_t& element, const word_t& word);
        bool sample(const element_t& element);
        void gibbs_sample(size_t n, size_t burnin);

        // constants
        ////////////////////////////////////////////////////////////////////////
        static const size_t BG_CLUSTER  = 0;
        static const size_t BG_LENGTH   = 1;
//        static const int TFBS_LENGTH = 42;
        static const size_t TFBS_LENGTH = 10;
        static const size_t NUCLEOTIDES = 4;

private:
        // data and clusters
        Data* data;
        ClusterManager* cluster_manager;

        // tags of special clusters
        cluster_tag_t bg_cluster_tag;

        // parameters
        double alpha;
        double lambda;

        // priors
        gsl_matrix* bg_alpha;
        gsl_matrix* tfbs_alpha;

        // gibbs sampler
        std::vector<std::vector<bool> > tfbs_start_positions;

        // gibbs sampler history
        size_t total_sampling_steps;
        std::vector<double> hist_switches;
        std::vector<double> hist_likelihood;
        std::vector<size_t> hist_num_clusters;
        std::vector<std::vector<double> > posterior;

        // keep track of the number of transcription factor binding sites
        size_t num_tfbs;
};

#endif /* DPM_HH */
