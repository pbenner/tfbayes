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
         DPM(double alpha, double lambda, size_t tfbs_length, const Data& data);
        ~DPM();

        // type definitions
        ////////////////////////////////////////////////////////////////////////
        typedef std::vector<std::vector<double> > posterior_t;
        typedef std::vector<std::vector<bool>   > tfbs_start_positions_t;

        // operators
        ////////////////////////////////////////////////////////////////////////
              Cluster& operator[](cluster_tag_t c)       { return (*_cluster_manager)[c]; }
        const Cluster& operator[](cluster_tag_t c) const { return (*_cluster_manager)[c]; }

        friend std::ostream& operator<<(std::ostream& o, const DPM& dpm);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t mixture_components() const;
        void mixture_weights(const word_t& word, double weights[], cluster_tag_t tags[]) const;
        const std::vector<std::vector<double> >& posterior() const {
                return _posterior;
        }
        const Data& data() const {
                return _data;
        }
        const ClusterManager& cluster_manager() const {
                return *_cluster_manager;
        }
        void add_word(const word_t& word, cluster_tag_t tag);
        void remove_word(const word_t& word, cluster_tag_t tag);
        void update_posterior(size_t sampling_steps);
        double compute_likelihood();
        bool valid_for_sampling(const element_t& element, const word_t& word);

        // constants
        ////////////////////////////////////////////////////////////////////////
        static const size_t BG_LENGTH   = 1;
               const size_t TFBS_LENGTH;
        static const size_t NUCLEOTIDES = 4;

private:
        // data and clusters
        const Data& _data;
        ClusterManager* _cluster_manager;

        // tags of special clusters
        cluster_tag_t bg_cluster_tag;

        // parameters
        double alpha;
        double lambda;

        // priors
        gsl_matrix* bg_alpha;
        gsl_matrix* tfbs_alpha;

        // record start positions of tfbs
        tfbs_start_positions_t tfbs_start_positions;

        // posterior distribution
        posterior_t _posterior;

        // keep track of the number of transcription factor binding sites
        size_t num_tfbs;
};

#endif /* DPM_HH */
