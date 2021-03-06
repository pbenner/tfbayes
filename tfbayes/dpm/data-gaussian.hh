/* Copyright (C) 2011-2013 Philipp Benner
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

#ifndef __TFBAYES_DPM_DATA_GAUSSIAN_HH__
#define __TFBAYES_DPM_DATA_GAUSSIAN_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <tfbayes/dpm/data.hh>
#include <tfbayes/dpm/datatypes.hh>
#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/dpm/indexer.hh>

class data_gaussian_t : public data_t<std::vector<double> >, public indexer_t {
        typedef data_t<std::vector<double> > base_t;
public:
        data_gaussian_t(size_t samples,
                        const std::matrix<double>& Sigma,
                        const std::vector<double>& pi);

        // operators
        ////////////////////////////////////////////////////////////////////////
        virtual data_gaussian_t& operator=(const data_i<std::vector<double> >& data) {
                using std::swap;
                data_gaussian_t tmp(static_cast<const data_gaussian_t&>(data));
                swap(*this, tmp);
                return *this;
        }
        // iterators
        ////////////////////////////////////////////////////////////////////////
        indexer_t::iterator begin() { return indices.begin(); }
        indexer_t::iterator end()   { return indices.end();   }

        indexer_t::const_iterator begin() const { return indices.begin(); }
        indexer_t::const_iterator end()   const { return indices.end();   }

        sampling_iterator sampling_begin() const
                { return sampling_indices.begin(); }
        sampling_iterator sampling_end() const
                { return sampling_indices.end();   }

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t elements() const;
        void shuffle();
        const std::matrix<double>& initial_means() const;
        const std::vector<double>& initial_cluster_assignments() const;

private:
        std::vector<index_t> indices;
        std::vector<index_t> sampling_indices;

        size_t _elements;
        size_t _cluster;

        // means for generating samples
        std::matrix<double> _mu;
        // initial cluster assignments
        std::vector<double> _initial_cluster_assignments;
};

#endif /* __TFBAYES_DPM_DATA_GAUSSIAN_HH__ */
