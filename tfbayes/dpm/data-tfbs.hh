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

#ifndef DATA_TFBS_HH
#define DATA_TFBS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <vector>
#include <string>

#include <gsl/gsl_matrix.h>

#include <code.hh>
#include <component-model.hh>
#include <indexer.hh>

class data_tfbs_t : public Indexer, public sequence_data_t<short> {
public:
         data_tfbs_t(const std::vector<std::string>& sequences, size_t tfbs_length);
         data_tfbs_t(const data_tfbs_t& data);
        ~data_tfbs_t();

        // iterators
        ////////////////////////////////////////////////////////////////////////

        // iterators over the full dataset (excluding masked nucleotides)
        iterator begin() { return indices.begin(); }
        iterator end()   { return indices.end();   }
        const_iterator begin() const { return indices.begin(); }
        const_iterator end()   const { return indices.end();   }

        // randomized iterator where in addition regions are excluded
        // where no binding site would fit
        sampling_iterator sampling_begin() const
                { return sampling_indices.begin(); }
        sampling_iterator sampling_end()   const
                { return sampling_indices.end();   }

        friend std::ostream& operator<< (std::ostream& o, const data_tfbs_t& data);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t elements() const { return _elements; };
        const std::vector<size_t>& sizes() const;
        void shuffle();

private:
        // all nucleotide positions in a vector (used for the gibbs sampler)
        std::vector<index_t*> indices;
        std::vector<index_t*> sampling_indices;

        // the raw nucleotide sequences
        const std::vector<std::string> sequences;
              std::vector<size_t     > sequences_length;

        // number of sequences and nucleotides
        size_t _n_sequences;
        size_t _elements;
};

#endif /* DATA_TFBS_HH */
