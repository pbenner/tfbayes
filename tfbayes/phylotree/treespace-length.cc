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
        // assert that we have two trees
        assert(tree_list.size() == 2);

        // convert trees
        ntree_t ntree1(tree_list.front());
        ntree_t ntree2(tree_list.back ());

        // destroy trees
        tree_list.front()->destroy();
        tree_list.back ()->destroy();

        // geodesic
        geodesic_t geodesic(ntree1, ntree2);

        cout << "length: "
             << geodesic.length()
             << endl;

        return 0;
}
