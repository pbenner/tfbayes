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

#include <algorithm>
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

static const size_t ALPHABET_SIZE = DpmTfbs::ALPHABET_SIZE; 
static const size_t PARSMM_DEPTH  = 2; 

// Standard output
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const nucleotide_sequence_t& sequence)
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
void add_subsequence(
        count_t* counts,
        size_t from, size_t to,
        const seq_context_t& context,
        const nucleotide_sequence_t& sequence)
{
        for(size_t i = from; i < min(sequence.size(), to+PARSMM_DEPTH+1); i++) {
                for(size_t c = 0; c < PARSMM_DEPTH+1; c++) {
                        if (context[i][c] != -1) {
                                counts[context[i][c]] += 1;
//                                cout << "Adding: " << i << ":" << context[i] << endl;
                        }
                }
        }
}

static
void remove_subsequence(
        count_t* counts,
        size_t from, size_t to,
        const seq_context_t& context,
        const nucleotide_sequence_t& sequence)
{
        for(size_t i = from; i < min(sequence.size(), to+PARSMM_DEPTH+1); i++) {
                for(size_t c = 0; c < PARSMM_DEPTH+1; c++) {
                        if (context[i][c] != -1) {
                                counts[context[i][c]] -= 1;
//                                cout << "Removing: " << i << ":" << context[i] << endl;
                        }
                }
        }
}

// Test functions
////////////////////////////////////////////////////////////////////////////////

static
void pt_init(static_pars_tree_t* pt) {
        for (size_t i = 0 ; i < pt->size * pt->as->size ; i++) {
                pt->dirichlet_params[i] = 1.0;
        }
}

static
double predictive(
        count_t* counts,
        size_t from, size_t to,
        static_pars_tree_t* pt,
        const seq_context_t& context,
        const nucleotide_sequence_t& sequence)
{
        double ml1 = pt_ln_marginal_likelihood(pt, counts);
        remove_subsequence(counts, from, to, context, sequence);
        double ml2 = pt_ln_marginal_likelihood(pt, counts);
        add_subsequence(counts, from, to, context, sequence);

        cout << "Marginal likelihood 1: P(y,x) = " << ml1 << endl;
        cout << "Marginal likelihood 2: P(x)   = " << ml2 << endl;

        return exp(ml1-ml2);
}

static void parsmm_test() __attribute__((unused));
static
void parsmm_test() {
        abstract_set_t* as = as_create(ALPHABET_SIZE);
        static_pars_tree_t* pt = pt_create(as, PARSMM_DEPTH);
        size_t counts_length = context_t::counts_size(ALPHABET_SIZE, PARSMM_DEPTH);;
        count_t* counts = (count_t*)malloc(counts_length*sizeof(count_t));

        /* init data structures */
        pt_init(pt);
        memset(counts, 0, counts_length*sizeof(count_t));

        /* test nucleotide sequence */
        const nucleotide_sequence_t sequence("GGGGGACGTCGATGCGTGATCGACTACGGCT");
        const seq_context_t context(sequence, PARSMM_DEPTH, ALPHABET_SIZE);
        cout << "Coded sequence:" << endl << sequence << endl << endl;
        cout << "Tree depth: "    << PARSMM_DEPTH << endl << endl;

        /**
         * test string:
         * GGGGGACGTCGATGCGTGATCGACTACGGCT
         *                |>---<|
         * remove:        GTGATCG
         * and obtain its probability from
         * the predictive distribution
         */
        add_subsequence(counts, 0, sequence.size(), context, sequence);
        print_counts(counts, counts_length);

        cout << "\nPredictive probability for subsequence GTGATCG at (15:21).\n"
             << "                        GGGGGACGTCGATGCGTGATCGACTACGGCT\n"
             << "                        0000013023012030201230132130032"
             << endl << endl;

        cout << "Predictive: "
             << predictive(counts, 15, 21, pt, context, sequence)
             << endl;

        (void)free(counts);
        pt_free(pt);
        as_free(as);
}

static void sanity_check() __attribute__((unused));
static void sanity_check() {
        abstract_set_t* as = as_create(ALPHABET_SIZE);
        static_pars_tree_t* pt = pt_create(as, PARSMM_DEPTH);
        size_t counts_length = context_t::counts_size(ALPHABET_SIZE, PARSMM_DEPTH);
        count_t* counts1 = (count_t*)malloc(counts_length*sizeof(count_t));
        count_t* counts2 = (count_t*)malloc(counts_length*sizeof(count_t));

        /* init data structures */
        pt_init(pt);
        memset(counts1, 0, counts_length*sizeof(count_t));
        memset(counts2, 0, counts_length*sizeof(count_t));

        const nucleotide_sequence_t sequence2("AGGGGACGTCGATGCGTGATCGACTACGGC");

        // use counts2 as reference
        const seq_context_t context2(sequence2, PARSMM_DEPTH, ALPHABET_SIZE);
        add_subsequence(counts2, 0, sequence2.size(), context2, sequence2);
        print_counts(counts2, counts_length);
        double ml_ref = pt_ln_marginal_likelihood(pt, counts2);
        double ml;

        // fill counts1 with sequences
        char n[4] = {'A', 'C', 'T', 'G'};
        string n_str("AGGGGACGTCGATGCGTGATCGACTACGGCXX");
        double sum = 0;
        for (short i = 0; i < (short)ALPHABET_SIZE; i++) {
                n_str[n_str.length()-2] = n[i];
                for (short j = 0; j < (short)ALPHABET_SIZE; j++) {
                        n_str[n_str.length()-1] = n[j];
                        const nucleotide_sequence_t sequence1(n_str);
                        const seq_context_t context1(sequence1, PARSMM_DEPTH, ALPHABET_SIZE);
                        add_subsequence(counts1, 0, sequence1.size(), context1, sequence1);

                        ml = pt_ln_marginal_likelihood(pt, counts1);
                        cout.precision(10);
                        cout << "Predictive: "
                             << "P(" << sequence1 << ")/P(" << sequence2 << ")"
                             << " = "
                             << ml << " - " << ml_ref << " = "
                             << exp(ml-ml_ref)
                             << endl;
                        sum += exp(ml-ml_ref);

                        remove_subsequence(counts1, 0, sequence1.size(), context1, sequence1);
                }
                cout << endl;
        }
        cout << "Sum: " << sum << endl;

        (void)free(counts1);
        (void)free(counts2);
        pt_free(pt);
        as_free(as);
}

int main()
{
        //parsmm_test();
        sanity_check();

        return 0;
}
