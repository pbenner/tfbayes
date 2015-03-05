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

const alphabet_t data_tfbs_t::alphabet = nucleotide_alphabet_t();

data_tfbs_t::data_tfbs_t(const string& phylogenetic_input)
        : sequence_data_t<code_t>(read_fasta(phylogenetic_input)),
          _n_sequences(size()),
          _elements(0)
{
        _complements = sequence_data_t<code_t>(*this);
        // compute complements
        for(size_t i = 0; i < _n_sequences; i++) {
                // loop over elements in a sequence
                for(size_t j = 0; j < operator[](i).size(); j++) {
                        // generate an index of this position
                        index_t index(i,j);
                        for (size_t k = 0; k < alphabet_size; k++) {
                                _complements[index][alphabet.complement(k)]
                                        = operator[](index)[k];
                        }
                }
        }
        // add indices
        for(size_t i = 0; i < _n_sequences; i++) {
                // loop over elements in a sequence
                for(size_t j = 0; j < operator[](i).size(); j++) {
                        // generate an index of this position
                        index_t index(i,j);
                        // push a new index to the list of indices
                        indices.push_back(index);
                        // push it also to the list of
                        // sampling indices
                        sampling_indices.push_back(index);
                }
                // increment the number of nucleotides
                _elements++;
        }
        shuffle();
}

void
data_tfbs_t::shuffle() {
        random_shuffle(sampling_indices.begin(), sampling_indices.end());
}

const sequence_data_t<__CODE_TYPE__>&
data_tfbs_t::complements() const
{
        return _complements;
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
                                 * larger than what is used for the
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
