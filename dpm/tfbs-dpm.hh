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

#ifndef TFBS_DPM_HH
#define TFBS_DPM_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include "data.hh"
#include "dpm.hh"
#include "cluster.hh"
#include "statistics.hh"
#include "tfbs-data.hh"

using namespace std;

class TfbsDPM : public DPM {
public:
        TfbsDPM(TfbsData* data);
        ~TfbsDPM();

        Distribution& posteriorPredictive(const Cluster::cluster& cluster);
        Distribution& predictive();

        double likelihood();
        void compute_statistics();

        TfbsData& get_data() {
                return *static_cast<TfbsData *>(da);
        }

        bool sample(Data::element& element);
        void gibbsSample(unsigned int steps);

        // constants
        static const int BG_LENGTH   = 1;
        static const int TFBS_LENGTH = 10;
        static const int NUCLEOTIDES = 4;

private:
        // parameters
        double lambda;

        // priors
        gsl_matrix* pd_tfbs_alpha;
        gsl_matrix* pd_bg_alpha;

        // distributions
        Distribution* predictiveDist_tfbs;
        Distribution* posteriorPredictiveDist_tfbs;
        Distribution* posteriorPredictiveDist_bg;

        // private methods
        void count_statistic(const Cluster::cluster& cluster, gsl_matrix* alpha, gsl_matrix* counts);
        bool check_element(Data::element& element);
        void assign_block(Data::element& element, Cluster::cluster_tag_t c);
        void release_block(Data::element& element);
};

#endif /* TFBS_DPM_HH */
