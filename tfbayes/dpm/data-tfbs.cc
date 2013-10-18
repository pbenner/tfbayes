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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <ctime>

#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include <tfbayes/fasta/fasta.hh>
#include <tfbayes/dpm/data-tfbs.hh>

using namespace std;

bool
data_tfbs_t::valid_sampling_index(const index_i& index, size_t tfbs_length) const
{
        // return false if the sequence ends before a tfbs could fit here
        if (operator[](index[0]).size() - index[1] < tfbs_length) {
                return false;
        }
        // check that all the elements in the range are not blank
        for (size_t i = 0; i < tfbs_length; i++) {
                seq_index_t tmp(index[0], index[1]+i);
                if (is_blank(tmp)) {
                        return false;
                }
        }
        // otherwise use this as a valid sampling index
        return true;
}

data_tfbs_t::data_tfbs_t(const string& phylogenetic_input, size_t tfbs_length)
        : sequence_data_t<code_t>(read_fasta(phylogenetic_input)),
          _n_sequences(size()),
          _elements(0)
{
        // loop over sequences
        for(size_t i = 0; i < _n_sequences; i++) {
                // loop over elements in a sequence
                for(size_t j = 0; j < operator[](i).size(); j++) {
                        // generate an index of this position
                        seq_index_t index(i,j);
                        // if there is a nucleotide at this position
                        if (!is_blank(index)) {
                                // push a new index to the list of indices
                                indices.push_back(new seq_index_t(index));
                                // if this index is also valid for
                                // sampling, i.e. a tfbs could fit here
                                if (valid_sampling_index(index, tfbs_length)) {
                                        // push it also to the list of
                                        // sampling indices
                                        sampling_indices.push_back(new seq_index_t(index));
                                }
                                // increment the number of nucleotides
                                _elements++;
                        }
                }
        }
        shuffle();
}

data_tfbs_t::data_tfbs_t(const data_tfbs_t& data)
        : sequence_data_t<code_t>(*this),
          _n_sequences(data._n_sequences),
          _elements(data._elements)
{
        for (indexer_t::const_iterator it = data.begin(); it != data.end(); it++) {
                index_i* index = (**it).clone();
                indices.push_back(index);
        }
        for (indexer_t::sampling_iterator it = data.sampling_begin();
             it != data.sampling_end(); it++)
        {
                index_i* index = (**it).clone();
                sampling_indices.push_back(index);
        }
        shuffle();
}

data_tfbs_t::~data_tfbs_t()
{
        for (indexer_t::iterator it = begin(); it != end(); it++) {
                delete(*it);
        }
        for (indexer_t::sampling_iterator it = sampling_begin();
             it != sampling_end(); it++) {
                delete(*it);
        }
}

void
data_tfbs_t::shuffle() {
        random_shuffle(sampling_indices.begin(), sampling_indices.end());
}

bool
data_tfbs_t::is_blank(const index_i& index) const
{
        // the position at index is blank if all counts of the
        // multinomial distribution are zero
        for (size_t i = 0; i < alphabet_size; i++) {
                if (operator[](index)[i] != 0.0) {
                        return false;
                }
        }
        return true;
}

#include <boost/regex.hpp>
#include <tfbayes/utility/strtools.hh>

sequence_data_t<data_tfbs_t::code_t>
data_tfbs_t::read_fasta(const string& file_name)
{
        sequence_data_t<data_tfbs_t::code_t> sequences;

        /* use a regular expression to match the multinomial counts of
         * a single column in the alignment, i.e. a line separated by
         * a semi-colon */
        boost::regex e("^"
                       "[[:space:]]*([[:digit:].]+)"
                       "[[:space:]]+([[:digit:].]+)"
                       "[[:space:]]+([[:digit:].]+)"
                       "[[:space:]]+([[:digit:].]+)"
                       "[[:space:]]+([[:digit:].]+)"
                       "[[:space:]]*"
                       "$");
        boost::smatch what;

        /* this automatically parses the fasta file format */
        FastaParser parser(file_name);

        /* storage for a single line */
        string line;

        /* the fasta parser returns a single line for each
         * sequence */
        while (parser) {
                line = parser();
                if (line == "") continue;
                /* store a single entry of the fasta file here */
                std::vector<data_tfbs_t::code_t> sequence;

                /* the count statistics for the multinomial
                 * distribution are separated by a semi-colon
                 */
                std::vector<std::string> result = token(line, ';');

                /* loop over all positions in a sequence */
                for (size_t i = 0; i < result.size(); i++) {
                        if(boost::regex_match(result[i], what, e, boost::match_extra))
                        {
                                /* the alphabet of the input can be
                                 * larger that what is used for the
                                 * sampler, for instance, gaps in the
                                 * alignment are used in the input as
                                 * fifth character, but are optional
                                 * for the sampler */
                                assert(what.size() > data_tfbs_t::alphabet_size);

                                data_tfbs_t::code_t entry;

                                for(size_t j = 0; j < data_tfbs_t::alphabet_size; j++)
                                {
                                        string submatch(what[j+1].first, what[j+1].second);
                                        entry[j] = atof(submatch.c_str());
                                }
                                sequence.push_back(entry);
                        }
                        else {
                                std::cerr << "WARNING: token `"
                                          << result[i]
                                          << "' did not match the regular expression."
                                          << std::endl;
                        }
                }
                /* save the result */
                sequences.push_back(sequence);
        }
        return sequences;
}

ostream& operator<< (ostream& o, const data_tfbs_t::code_t& entry)
{
        for (size_t k = 0; k < data_tfbs_t::alphabet_size; k++)
        {
                o.precision(8);
                o.width(10);
                o << entry[k]
                  << " ";
        }
        o << ";" << endl;

        return o;
}

/* Print this object in a fasta file format */
ostream& operator<< (ostream& o, const data_tfbs_t& data)
{
        for (size_t i = 0; i < data.size(); i++)
        {
                o << ">sequence "<< i << endl;

                for (size_t j = 0; j < data[i].size(); j++)
                {
                        o << data[i][j];
                }
        }

        return o;
}
