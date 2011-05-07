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

#include <tfbayes/logarithmetic.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#include "cluster.hh"
#include "data.hh"
#include "dpm.hh"
#include "statistics.hh"

using namespace std;

DPM::DPM(Data* data) : da(data), cl(*da), alpha(1.0) {
        hist_switches.push_back(0);
        hist_likelihood.push_back(likelihood());
}

DPM::~DPM() {
        delete(da);
}

bool DPM::sample(Data::element& element) {
        Cluster::cluster_tag_t old_cluster_tag = cl.getClusterTag(element);
        cl.release(element);
        Distribution& pred = predictive();
        Cluster::size_type num_clusters = cl.size();
        double ln_weights[num_clusters+1];
        double weights[num_clusters+1];
        Cluster::cluster_tag_t tags[num_clusters+1];
        double log_sum = -HUGE_VAL;

        // compute weights of existing clusters
        Cluster::cluster_tag_t i = 0;
        for (Cluster::iterator it = cl.begin(); it != cl.end(); it++) {
                Distribution& postPred = posteriorPredictive(**it);
                double num_elements    = (double)(*it)->elements.size();
                ln_weights[i] = log(num_elements) + postPred.ln_pdf(element.x);
                tags[i]       = (*it)->tag;
                // normalization constant
                log_sum = logadd(ln_weights[i], log_sum);
                i++;
        }
        // add the tag of a new class and compute their weight
        printf("alpha: %f\n", alpha);
        ln_weights[num_clusters] = log(alpha) + pred.ln_pdf(element.x);
        tags[num_clusters]       = cl.next_free_cluster()->tag;
        log_sum = logadd(ln_weights[num_clusters], log_sum);

        // normalize
        for (i = 0; i < num_clusters+1; i++) {
                printf("ln_weights[%ld]: %f\n", i, ln_weights[i]);
                printf("log_sum: %f\n", log_sum);
                weights[i] = expl(ln_weights[i] - log_sum);
                printf("weights[%ld]: %f\n", i, weights[i]);
        }

        // draw a new cluster for the element
        gsl_ran_discrete_t* gdd = gsl_ran_discrete_preproc(num_clusters+1, weights);
        i = gsl_ran_discrete(_r, gdd);
        gsl_ran_discrete_free(gdd);

        cl.assign(element, tags[i]);
        printf("\n\n");

        return old_cluster_tag != tags[i];
}

void DPM::gibbsSample(unsigned int steps) {
        for (unsigned int i = 0; i < steps; i++) {
                double sum = 0;
                for (Data::iterator it = da->begin(); it != da->end(); it++) {
                        bool switched = sample(*it);
                        if (switched) sum+=1;
                }
                hist_switches.push_back(sum/da->size());
                compute_statistics();
        }
}

ostream& operator<< (ostream& o, DPM const& dpm)
{
        return o << dpm.cl;
}
