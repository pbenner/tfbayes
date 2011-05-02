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

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "tfbs-dpm.hh"

using namespace std;

TfbsDPM::TfbsDPM(TfbsData* data)
        : DPM(data)
{
        // initialize distributions
//        predictiveDist          = new BivariateNormal(predictive_cov, mu_0);
//        posteriorPredictiveDist = new BivariateNormal();
}

TfbsDPM::~TfbsDPM() {
        delete(predictiveDist);
        delete(posteriorPredictiveDist);
}

Distribution& TfbsDPM::posteriorPredictive(const Cluster::cluster& cluster) {

        return *posteriorPredictiveDist;
}

Distribution& TfbsDPM::predictive() {
        return *predictiveDist;
}

double TfbsDPM::likelihood() {
        return 0.0;
}

void TfbsDPM::compute_statistics() {
}
