/* Copyright (C) 2012 Philipp Benner
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

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>

#include <tfbayes/phylotree/interface.hh>
#include <tfbayes/phylotree/utility.hh>

using namespace std;

__BEGIN_DECLS

// library interface
////////////////////////////////////////////////////////////////////////////////

alignment_t<code_t>* alignment_new(size_t length, pt_root_t* pt_root)
{
        // initialize the alignment with -1, meaning that positions
        // that are not later initialized otherwise will be ignored by
        // the phylogenetic model; i.e. if a species is not part of
        // the multiple alignment it will be fully ignored
        return new alignment_t<code_t>(length, -1, *pt_root);
}

void alignment_set(alignment_t<code_t>* alignment, const char* taxon, vector_t* record)
{
        assert(alignment->length() == record->size);

        for (size_t i = 0; i < alignment->length(); i++) {
                alignment->operator[](taxon)[i] = (code_t)record->vec[i];
        }
}

vector_t* alignment_marginal_likelihood(const pt_root_t* tree, alignment_t<code_t>* alignment, vector_t* prior)
{
        exponent_t<code_t, alphabet_size> alpha;
        vector_t* result = alloc_vector(alignment->length());

        for (size_t i = 0; i < alphabet_size; i++) {
                alpha[i] = prior->vec[i];
        }

        /* go through the alignment and compute the marginal
         * likelihood for each position */
        for (alignment_t<code_t>::iterator it = alignment->begin(); it != alignment->end(); it++) {
                result->vec[it.position()] = pt_marginal_likelihood<code_t, alphabet_size>(*tree, *it, alpha);
        }
        return result;
}

vector_t* alignment_scan(const pt_root_t* tree, alignment_t<code_t>* alignment, matrix_t* c_counts)
{
        vector_t* result = alloc_vector(alignment->length());
        vector<exponent_t<code_t, alphabet_size> > counts;

        for (size_t j = 0; j < c_counts->columns; j++) {
                exponent_t<code_t, alphabet_size> tmp;
                for (size_t i = 0; i < alphabet_size; i++) {
                        // check that all counts are positive
                        assert(c_counts->mat[i][j] >= 0.0);
                        tmp[i] = c_counts->mat[i][j];
                }
                counts.push_back(tmp);
        }

        for (alignment_t<code_t>::iterator it = alignment->begin(); it != alignment->end(); it++) {
                result->vec[it.position()] = 0;
                if (it.position() + counts.size() > alignment->length()) {
                        // do not exit the loop here so that every
                        // position is initialized
                        continue;
                }
                for (alignment_t<code_t>::iterator is(it); is.position() < it.position() + counts.size(); is++) {
                        result->vec[it.position()] += pt_marginal_likelihood<code_t, alphabet_size>(*tree, *is, counts[is.position()-it.position()]);
                }
        }
        return result;
}

void alignment_free(alignment_t<code_t>* alignment)
{
        delete(alignment);
}

__END_DECLS
