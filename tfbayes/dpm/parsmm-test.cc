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

#include <sstream>
#include <string.h>
#include <math.h>

#include <abysmal-stack.hh>
#include <dpm-tfbs.hh>
#include <nucleotide-sequence.hh>

#include <tfbayes/logarithmetic.h>
#include <tfbayes/fastlog.h>

#include <parsmm/abstract_set.h>
#include <parsmm/static_pars_tree.h>

using namespace std;

// Constants
////////////////////////////////////////////////////////////////////////////////

static const size_t ALPHABET_SIZE = DPM_TFBS::ALPHABET_SIZE; 
static const size_t PARSMM_DEPTH  = DPM_TFBS::PARSMM_DEPTH; 

// Standard output
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const AbysmalStack<count_t>& stack)
{
        for (size_t i = 0; i < PARSMM_DEPTH+1; i++) {
                o << stack[i];
        }

        return o;
}

ostream& operator<< (ostream& o, const NucleotideSequence& sequence)
{
        for (size_t i = 0; i < sequence.size(); i++) {
                o << sequence[i];
        }

        return o;
}

static
void print_counts(count_t* counts, size_t length)
{
        cout << "Counts:" << endl;
        for (size_t i = 0; i < length; i++) {
                cout << counts[i] << " ";
        }
        cout << endl;
}

// Counts
////////////////////////////////////////////////////////////////////////////////

static
size_t counts_position(const AbysmalStack<count_t>& stack)
{
        size_t position = pow(ALPHABET_SIZE, PARSMM_DEPTH+1) - 1;

        for (size_t i = 0; i < PARSMM_DEPTH+1; i++) {
                position -= stack[i]*pow(ALPHABET_SIZE, i);
        }

        return position;
}

static
void compute_counts(count_t* counts, const NucleotideSequence& sequence)
{
        AbysmalStack<count_t> stack(PARSMM_DEPTH+1);

        cout << "Counts positions:" << endl;
        for (size_t i = 0; i < sequence.size(); i++) {
                stack.push(sequence[i]);
                if (i >= PARSMM_DEPTH) {
                        size_t pos = counts_position(stack);
                        counts[pos] += 1;
                        cout << stack << " -> " << pos << endl;
                }
        }
        cout << endl;
}

// Test functions
////////////////////////////////////////////////////////////////////////////////

static
void parsmm_test() {
        abstract_set_t* as = as_create(ALPHABET_SIZE);
        static_pars_tree_t* pt = pt_create(as, PARSMM_DEPTH);
        AbysmalStack<count_t> stack(PARSMM_DEPTH+1);
        size_t counts_length = pow(ALPHABET_SIZE, PARSMM_DEPTH+1);
        count_t* counts = (count_t*)malloc(counts_length*sizeof(count_t));
        memset(counts, 0, counts_length*sizeof(count_t));

        /* test nucleotide sequence */
        const NucleotideSequence sequence("GGGGGACGTCGATGCGTGATCGACTACGGCT");
        cout << "Coded sequence:" << endl << sequence << endl << endl;

        cout << "Tree depth: " << PARSMM_DEPTH << endl << endl;

        /**
         * test string:
         * GGGGGACGTCGATGCGTGATCGACTACGGCT
         *                |>---<|
         * remove:        GTGATCG
         * and obtain its probability from
         * the predictive distribution
         */
        compute_counts(counts, sequence);
        print_counts(counts, counts_length);

        cout << "Marginal likelihood: "
             << pt_ln_marginal_likelihood(pt, counts)
             << endl;

        cout << "Removing subsequence GTGATCG at (15:21)."
             << endl << endl;

        (void)free(counts);
}

int main()
{
        parsmm_test();

        return 0;
}
