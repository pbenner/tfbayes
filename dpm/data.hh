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

#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>

#include <datatypes.hh>

class Data {
public:
        Data(size_t n, char *sequences[], cluster_tag_t default_tag);
        ~Data();

        // type definitions
        typedef std::vector<element_t >::iterator iterator;
        typedef std::vector<element_t >::const_iterator const_iterator;

        typedef std::vector<element_t*>::iterator iterator_randomized;
        typedef std::vector<element_t*>::const_iterator const_iterator_randomized;

        // iterators
        ////////////////////////////////////////////////////////////////////////
        iterator begin() { return elements.begin(); }
        iterator end()   { return elements.end();   }

        const_iterator begin() const { return elements.begin(); }
        const_iterator end()   const { return elements.end(); }

        iterator_randomized begin_randomized() { return elements_randomized.begin(); }
        iterator_randomized end_randomized()   { return elements_randomized.end(); }

        const_iterator_randomized begin_randomized() const
                { return elements_randomized.begin(); }
        const_iterator_randomized end_randomized() const
                { return elements_randomized.end(); }

        // operators
        ////////////////////////////////////////////////////////////////////////
              element_t& operator[](size_t i);
        const element_t& operator[](size_t i) const;

//        friend ostream& operator<< (std::ostream& o, Data const& data);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t size() { return elements.size(); }
        void shuffle();

        int  num_successors(const element_t& e) const;

        bool valid_for_sampling(const element_t& element, size_t length, word_t& word);

        size_t get_n_sequences() {
                return n_sequences;
        }
        size_t get_sequence_length(size_t i) {
                return sequences_length[i];
        }

        void record_cluster_assignment(const word_t& word, cluster_tag_t tag) {
                for (size_t i = 0; i < word.length; i++) {
                        cluster_assignments[word.sequence][word.position+i] = tag;
                }
        }

        cluster_tag_t get_cluster_tag(const element_t& element) const {
                return cluster_assignments[element.sequence][element.position];
        }

private:
        // all nucleotide positions in a vector (used for the gibbs sampler)
        std::vector<element_t > elements;
        std::vector<element_t*> elements_randomized;

        // the raw nucleotide sequences
        std::vector<std::string> sequences;
        std::vector<size_t> sequences_length;
        size_t n_sequences;

        // assignments to clusters
        std::vector<std::vector<cluster_tag_t> > cluster_assignments;
};

//ostream& operator<< (ostream& o, element_t const& element);

#endif /* DATA_HH */
