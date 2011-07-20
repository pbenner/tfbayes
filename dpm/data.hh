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

#ifndef DATA_HH
#define DATA_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <vector>
#include <string>

#include <gsl/gsl_matrix.h>

#include <code.hh>
#include <datatypes.hh>

class Data_TFBS : public sequence_data_t<short> {
public:
         Data_TFBS(const std::vector<std::string>& sequences);

        // type definitions
        typedef std::vector<index_t>::iterator iterator;
        typedef std::vector<index_t>::const_iterator const_iterator;

        typedef std::vector<index_t*>::iterator iterator_randomized;
        typedef std::vector<index_t*>::const_iterator const_iterator_randomized;

        // iterators
        ////////////////////////////////////////////////////////////////////////
        iterator begin() { return indices.begin(); }
        iterator end()   { return indices.end(); }

        const_iterator begin() const { return indices.begin(); }
        const_iterator end()   const { return indices.end(); }

        const_iterator_randomized begin_randomized() const
                { return indices_randomized.begin(); }
        const_iterator_randomized end_randomized() const
                { return indices_randomized.end(); }

        // operators
        ////////////////////////////////////////////////////////////////////////
        const index_t& operator[](size_t i) const;

        friend std::ostream& operator<< (std::ostream& o, const Data_TFBS& data);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t size() const;
        size_t length() const;
        size_t length(size_t i) const;
        const std::vector<size_t>& lengths() const;
        void shuffle();

private:
        // all nucleotide positions in a vector (used for the gibbs sampler)
        std::vector<index_t > indices;
        std::vector<index_t*> indices_randomized;

        // the raw nucleotide sequences
        const std::vector<std::string> sequences;
              std::vector<size_t     > sequences_length;

        // number of sequences and nucleotides
        size_t _n_sequences;
        size_t _size;
};

//std::ostream& operator<< (std::ostream& o, const word_t& word);

#endif /* DATA_HH */
