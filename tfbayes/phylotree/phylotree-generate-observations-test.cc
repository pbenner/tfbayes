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

#include <iostream>
#include <iomanip>
#include <vector>
#include <sys/time.h>

#include <boost/format.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/phylotree-generate-observations.hh>
#include <tfbayes/uipac/code.hh>

using namespace std;
using boost::format;

#define alphabet_size 5
typedef short code_t;

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

void print_alignment(const pt_root_t* root, const alignment_t<code_t>& alignment)
{
        for (size_t i = 0; i < alignment.n_species(); i++) {
                cout << format("%20s: ")
                        % root->leaves[i]->name;

                for (size_t j = 0; j < alignment.length(); j++) {
                        cout << decode_nucleotide(alignment[i][j]);
                }
                cout << endl;
        }
}

int main(void)
{
        // number of columns in the alignment
        const size_t N = 100;

        init();

        // gsl random number generator
        gsl_rng_env_setup();
        const gsl_rng_type * T = gsl_rng_default;
        gsl_rng * r = gsl_rng_alloc(T);

        // parse tree
        list<pt_root_t*> tree_list = parse_tree_list();
        assert(tree_list.size() == 1);
        pt_root_t* pt_root = tree_list.front();

        // alignment
        alignment_t<code_t> alignment(N, -1, pt_root);

        // pseudo counts
        vector<double> alpha(alphabet_size, 0);
        double a = 10.0;
        alpha[0] = a;
        alpha[1] = a;
        alpha[2] = a;
        alpha[3] = a;
        alpha[4] = 20;

        // generate
        for (size_t i = 0; i < N; i++) {
                vector<double> stationary   = dirichlet_sample<alphabet_size>(alpha, r);
                vector<code_t> observations = pt_generate_observations<code_t, alphabet_size>(pt_root, stationary);

                for (size_t k = 0; k < observations.size(); k++) {
                        alignment[k][i] = observations[k];
                }
        }
        // print result
        print_alignment(pt_root, alignment);

        // clean up
        pt_root->destroy();
        gsl_rng_free (r);
}
