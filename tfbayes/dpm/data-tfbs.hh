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

         data_tfbs_t(const std::string& phylogenetic_input, size_t tfbs_length);
         data_tfbs_t(const data_tfbs_t& data);
        ~data_tfbs_t();

        friend void swap(data_tfbs_t& first, data_tfbs_t& second) {
                using std::swap;
                swap(static_cast<sequence_data_t<__CODE_TYPE__>&>(first),
                     static_cast<sequence_data_t<__CODE_TYPE__>&>(second));
                swap(first.indices,          second.indices);
                swap(first.sampling_indices, second.sampling_indices);
                swap(first._n_sequences,     second._n_sequences);
                swap(first._elements,        second._elements);
        }

        // operators
        ////////////////////////////////////////////////////////////////////////
        virtual data_tfbs_t& operator=(const data_i<code_t>& data) {
                using std::swap;
                data_tfbs_t tmp(static_cast<const data_tfbs_t&>(data));
                swap(*this, tmp);
                return *this;
        }
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
        // true if there is no nucleotide at position given by index
        bool is_blank(const index_i& index) const;
        size_t elements() const { return _elements; };
        void shuffle();
        bool valid_sampling_index(const index_i& index, size_t tfbs_length) const;
        const sequence_data_t<__CODE_TYPE__>& complements() const;

        static sequence_data_t<data_tfbs_t::code_t> read_fasta(const std::string& file_name);

private:
        // all nucleotide positions in a vector (used for the gibbs sampler)
        std::vector<index_i*> indices;
        std::vector<index_i*> sampling_indices;
        // complements
        sequence_data_t<__CODE_TYPE__> _complements;

        // number of sequences and nucleotides
        size_t _n_sequences;
        size_t _elements;
};

std::ostream& operator<< (std::ostream& o, const data_tfbs_t::code_t& entry);

#endif /* __TFBAYES_DPM_DATA_TFBS_HH__ */
