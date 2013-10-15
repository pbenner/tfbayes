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

#include <tfbayes/phylotree/phylotree-parsetree.hh>
#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/alignment/sequence.hh>

using namespace std;

int
main(void)
{
        pt_root_t tree = parse_tree_list("alignment-test.nh").front();
        alignment_t<> alignment("alignment-test.fa", tree);

        cout << print_alignment_pretty(alignment) << endl;
        cout << print_alignment_fasta (alignment) << endl;

        alignment_set_t<> alignment_set("alignment-test.mfa", tree);

        for (alignment_set_t<>::const_iterator it = alignment_set.begin();
             it != alignment_set.end(); it++) {
                cout << print_alignment_pretty(*it) << endl;
        }

        return 0;
}
