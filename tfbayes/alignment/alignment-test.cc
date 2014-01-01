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

#include <iostream>

#include <tfbayes/phylotree/parsetree.hh>
#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/alignment/sequence.hh>

using namespace std;

int
main(void)
{
        // alignment with phylogenetic tree
        pt_root_t tree = parse_tree_list("alignment-test.nh").front();
        alignment_t<> alignment1("alignment-test.fa", tree);

        cout << print_alignment_pretty(alignment1) << endl;
        cout << print_alignment_fasta (alignment1) << endl;

        alignment_set_t<> alignment_set1("alignment-test.mfa", tree);

        for (alignment_set_t<>::const_iterator it = alignment_set1.begin();
             it != alignment_set1.end(); it++) {
                cout << print_alignment_pretty(*it) << endl;
        }

        // alignment with phylogenetic tree
        alignment_t<> alignment2("alignment-test.fa");
        cout << print_alignment_pretty(alignment2) << endl;

        alignment_set_t<> alignment_set2("alignment-test.mfa");
        for (alignment_set_t<>::const_iterator it = alignment_set2.begin();
             it != alignment_set2.end(); it++) {
                cout << print_alignment_pretty(*it) << endl;
        }

        return 0;
}
