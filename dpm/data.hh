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

#include "data.hh"

using namespace std;

class Data {
public:
        Data(int n, int m, char *sequences[], int *clusters[]);
        ~Data();

        // type definitions
        typedef size_t type;
        typedef vector<type> x_t;
        typedef struct {
                x_t x;                // value of the element
                int tag;              // tag to identify the element
                                      // (used for cluster assignments)
                int original_cluster; // from which cluster this
                                      // elemenent was originally sampled
        } element;
        typedef vector<element>::size_type size_type;

        typedef vector<element >::iterator iterator;
        typedef vector<element >::const_iterator const_iterator;

        typedef vector<element*>::iterator iterator_randomized;
        typedef vector<element*>::const_iterator const_iterator_randomized;

        // iterators
        ////////////////////////////////////////////////////////////////////////
        iterator begin() { return elements.begin(); }
        iterator end()   { return elements.end(); }

        const_iterator begin() const { return elements.begin(); }
        const_iterator end()   const { return elements.end(); }

        iterator_randomized begin_randomized() { return elements_randomized.begin(); }
        iterator_randomized end_randomized()   { return elements_randomized.end(); }

        const_iterator_randomized begin_randomized() const
                { return elements_randomized.begin(); }
        const_iterator_randomized end_randomized()   const
                { return elements_randomized.end(); }


        // operators
        ////////////////////////////////////////////////////////////////////////
              element& operator[](size_type i);
        const element& operator[](size_type i) const;

        friend ostream& operator<< (std::ostream& o, Data const& data);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_type size() { return elements.size(); }
        void shuffle();

        iterator find(const element& elem);
        char get_nucleotide(const element& e) const;
        void get_nucleotide(const element& e, int n, char *buf) const;
        int num_successors(const element& e);

private:
        vector<element > elements;
        vector<element*> elements_randomized;

        vector<string> sequences;
        size_t n_sequences;
        size_t sequence_length;
};

ostream& operator<< (ostream& o, Data::element const& element);

#endif /* DATA_HH */
