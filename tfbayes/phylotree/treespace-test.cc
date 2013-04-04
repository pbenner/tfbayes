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

#include <treespace.hh>

using namespace std;

void test1()
{
        // prepare splits
        set<size_t> s1;
        s1.insert(3);
        s1.insert(4);
        set<size_t> s2;
        s2.insert(1);
        s2.insert(2);
        set<size_t> s3;
        s3.insert(0);
        s3.insert(5);
        nsplit_t e1(5, s1);
        nsplit_t e2(5, s2);
        nsplit_t e3(5, s3);
        // construct tree
        vector<nsplit_t> splits;
        vector<double>   int_d;
        // internal edge lengths
        splits.push_back(e1); int_d.push_back(5.5);
        splits.push_back(e2); int_d.push_back(2.5);
        splits.push_back(e3); int_d.push_back(1.5);
        // leaf edge lengths
        vector<double> leaf_d(6, 1);
        leaf_d[0] = 1.5;
        leaf_d[1] = 2.5;
        leaf_d[2] = 3.5;
        leaf_d[3] = 4.5;
        leaf_d[4] = 5.5;
        // leaf names
        vector<string> leaf_names;
        leaf_names.push_back("leaf 1");
        leaf_names.push_back("leaf 2");
        leaf_names.push_back("leaf 3");
        leaf_names.push_back("leaf 4");
        leaf_names.push_back("leaf 5");
        leaf_names.push_back("leaf 6");
        ntree_t ntree(splits, int_d, leaf_d, leaf_names);
        // print tree
        pt_root_t* tree = ntree.export_tree();
        cout << "exporting tree:" << endl;
        tree->print_phylotree(cout);
        tree->destroy();
}

int
main(void)
{
        test1();

        return 0;
}
