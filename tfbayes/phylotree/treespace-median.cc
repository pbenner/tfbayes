/* Copyright (C) 2013 Philipp Benner
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

#include <cstdio>

#include <tfbayes/phylotree/treespace.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>

using namespace std;

int
main(void)
{
        // parse trees from file
        list<pt_root_t*> tree_list = parse_tree_list();
        list<ntree_t   > ntree_list;
        // convert trees
        for (list<pt_root_t*>::const_iterator it = tree_list.begin();
             it != tree_list.end(); it++) {
                pt_root_t* pt_root = *it;

                ntree_list.push_back(pt_root);
                pt_root->destroy();
        }

        cout << "Computing median" << endl
             << "------------------------------------------------------------"
             << endl;

        cout << median_tree(ntree_list)
             << endl;

        return 0;
}
