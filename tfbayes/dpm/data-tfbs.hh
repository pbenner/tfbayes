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

#ifndef __TFBAYES_DPM_DATA_TFBS_HH__
#define __TFBAYES_DPM_DATA_TFBS_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <vector>
#include <string>

#include <boost/array.hpp>

#include <tfbayes/dpm/data.hh>
#include <tfbayes/dpm/indexer.hh>

#include <tfbayes/uipac/alphabet.hh>

#define __ALPHABET_SIZE__ 5
#define __CODE_TYPE__ boost::array<double, __ALPHABET_SIZE__>

/* data_tfbs_t is the container for the data that provides some basic
 * operations like indexing and iteration */
class data_tfbs_t : public indexer_t, public sequence_data_t<__CODE_TYPE__> {
public:
        // typedefs
        ////////////////////////////////////////////////////////////////////////

        /* This type defines a singe entry in a sequence. Since the basic
         * type is an approximate multinomial distribution, we need to store
         * the count statistics for each nucleotide. */
        static const size_t alphabet_size = __ALPHABET_SIZE__;
        static const alphabet_t alphabet;
        typedef __CODE_TYPE__ code_t;

        // constructors
        ////////////////////////////////////////////////////////////////////////

        data_tfbs_t(const std::string& phylogenetic_input);

        // iterators
        ////////////////////////////////////////////////////////////////////////

        // iterators over the full dataset (excluding masked nucleotides)
        indexer_t::iterator begin() { return indices.begin(); }
        indexer_t::iterator end()   { return indices.end();   }
        indexer_t::const_iterator begin() const { return indices.begin(); }
        indexer_t::const_iterator end()   const { return indices.end();   }

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
        void shuffle();
        const sequence_data_t<__CODE_TYPE__>& complements() const;

        static sequence_data_t<data_tfbs_t::code_t> read_fasta(const std::string& file_name);

private:
        // all nucleotide positions in a vector (used for the gibbs sampler)
        std::vector<index_t> indices;
        std::vector<index_t> sampling_indices;
        // complements
        sequence_data_t<__CODE_TYPE__> _complements;

        // number of sequences and nucleotides
        size_t _n_sequences;
        size_t _elements;
};

std::ostream& operator<< (std::ostream& o, const data_tfbs_t::code_t& entry);

#endif /* __TFBAYES_DPM_DATA_TFBS_HH__ */
