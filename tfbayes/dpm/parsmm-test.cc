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

// Context
////////////////////////////////////////////////////////////////////////////////

class Context : public std::vector<ssize_t>
{
public:
        Context(const nucleotide_sequence_t& sequence, size_t depth, size_t alphabet_size) {
                AbysmalStack<count_t> stack(depth+1);
                size_t position;

                for (size_t i = 0; i < sequence.size(); i++) {
                        stack.push(sequence[i]);
                        if (i >= depth) {
                                position = pow(alphabet_size, depth+1) - 1;
                                for (size_t i = 0; i < depth+1; i++) {
                                        position -= stack[i]*pow(alphabet_size, i);
                                }
                                cout << stack << " -> " << position << endl;
                                push_back(position);
                        }
                        else {
                                push_back(-1);
                        }
                }
        }
};

// Counts
////////////////////////////////////////////////////////////////////////////////

static
void add_subsequence(
        count_t* counts,
        size_t from, size_t to,
        const Context& context,
        const nucleotide_sequence_t& sequence)
{
        for(size_t i = from; i < min(sequence.size(), to+PARSMM_DEPTH+1); i++) {
                if (context[i] != -1) {
                        counts[context[i]] += 1;
                        cout << "Adding: " << i << ":" << context[i] << endl;
                }
        }
}

static
void remove_subsequence(
        count_t* counts,
        size_t from, size_t to,
        const Context& context,
        const nucleotide_sequence_t& sequence)
{
        for(size_t i = from; i < min(sequence.size(), to+PARSMM_DEPTH+1); i++) {
                if (context[i] != -1) {
                        counts[context[i]] -= 1;
                        cout << "Removing: " << i << ":" << context[i] << endl;
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
        const Context& context,
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

static
void parsmm_test() {
        abstract_set_t* as = as_create(ALPHABET_SIZE);
        static_pars_tree_t* pt = pt_create(as, PARSMM_DEPTH);
        size_t counts_length = pow(ALPHABET_SIZE, PARSMM_DEPTH+1);
        count_t* counts = (count_t*)malloc(counts_length*sizeof(count_t));

        /* init data structures */
        pt_init(pt);
        memset(counts, 0, counts_length*sizeof(count_t));

        /* test nucleotide sequence */
        const nucleotide_sequence_t sequence("GGGGGACGTCGATGCGTGATCGACTACGGCT");
//        const nucleotide_sequence_t sequence("AAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
        const Context context(sequence, PARSMM_DEPTH, ALPHABET_SIZE);
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

int main()
{
        parsmm_test();

        return 0;
}
