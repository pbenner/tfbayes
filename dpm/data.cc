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

#include <iostream>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <ctime>

#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include <clustermanager.hh>
#include <data.hh>
#include <distribution.hh>

using namespace std;

Data::Data(size_t n, char *sequences[])
        : n_sequences(n)
{
        for(size_t i = 0; i < n; i++) {
                size_t m = strlen(sequences[i]);
                this->sequences.push_back(sequences[i]);
                this->sequences_length.push_back(m);
                for(size_t j = 0; j < m; j++) {
                        element_t e = {i, j};
                        elements.push_back(e);
                }
        }
        for (Data::iterator it = begin(); it != end(); it++) {
                elements_randomized.push_back(&(*it));
        }
        shuffle();
}

Data::~Data() {
}

void
Data::shuffle() {
        random_shuffle(elements_randomized.begin(), elements_randomized.end());
}

const element_t&
Data::operator[](size_t i) const {
        return elements[i];
}

element_t&
Data::operator[](size_t i) {
        return elements[i];
}

char
Data::operator[](element_t element) const {
        return sequences[element.sequence][element.position];
}

const word_t
Data::get_word(const element_t& element, size_t length) const {
        const size_t sequence = element.sequence;
        const size_t position = element.position;
        word_t word = {sequence, position, length, sequences};

        return word;
}

size_t
Data::length() const {
        return n_sequences;
}

size_t
Data::length(size_t i) const {
        if (i < n_sequences) {
                return sequences_length[i];
        }
        else {
                return 0;
        }
}

ostream&
operator<< (ostream& o, const word_t& word) {
        o << "(" << word.sequence << ":" << word.position << ":";
        for (size_t i = 0; i < word.length; i++) {
                o << word.sequences[word.sequence][word.position+i];
        }
        o << ")";

        return o;
}

ostream& operator<< (ostream& o, const Data& data)
{
        for (size_t i = 0; i < data.length(); i++) {
                o << data.sequences[i] << endl;
        }
        return o;
}
