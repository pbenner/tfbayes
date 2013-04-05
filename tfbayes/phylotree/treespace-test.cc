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
        nedge_t e1(5, s1, 5.5);
        nedge_t e2(5, s2, 2.5);
        nedge_t e3(5, s3, 1.5);
        // construct tree
        nedge_set_t nedge_set;
        // internal edge lengths
        nedge_set.push_back(e1);
        nedge_set.push_back(e2);
        nedge_set.push_back(e3);
        // leaf edge lengths
        vector<double> leaf_d(6, 1);
        leaf_d[0] = 0.5;
        leaf_d[1] = 1.5;
        leaf_d[2] = 2.5;
        leaf_d[3] = 3.5;
        leaf_d[4] = 4.5;
        leaf_d[5] = 5.5;
        // leaf names
        vector<string> leaf_names;
        leaf_names.push_back("leaf 0");
        leaf_names.push_back("leaf 1");
        leaf_names.push_back("leaf 2");
        leaf_names.push_back("leaf 3");
        leaf_names.push_back("leaf 4");
        leaf_names.push_back("leaf 5");
        ntree_t ntree(nedge_set, leaf_d, leaf_names);
        // print tree
        pt_root_t* tree = ntree.export_tree();
        cout << "exporting tree:" << endl;
        tree->print_phylotree(cout);
        tree->destroy();
}

void test2()
{
        // prepare splits
        set<size_t> s1;
        s1.insert(3);
        s1.insert(4);
        s1.insert(5);
        set<size_t> s2;
        s2.insert(4);
        s2.insert(5);
        set<size_t> s3;
        s3.insert(6);
        s3.insert(7);
        set<size_t> s4;
        s4.insert(0);
        s4.insert(1);
        s4.insert(2);
        set<size_t> s5;
        s5.insert(0);
        s5.insert(1);
        nedge_t e1(7, s1, 1);
        nedge_t e2(7, s2, 1);
        nedge_t e3(7, s3, 1);
        nedge_t e4(7, s4, 1);
        nedge_t e5(7, s5, 1);
        // construct tree
        nedge_set_t nedge_set;
        // internal edge lengths
        nedge_set.push_back(e1);
        nedge_set.push_back(e2);
        nedge_set.push_back(e3);
        nedge_set.push_back(e4);
        nedge_set.push_back(e5);
        // leaf edge lengths
        vector<double> leaf_d(8, 1);
        leaf_d[0] = 0.5;
        leaf_d[1] = 1.5;
        leaf_d[2] = 2.5;
        leaf_d[3] = 3.5;
        leaf_d[4] = 4.5;
        leaf_d[5] = 5.5;
        leaf_d[6] = 6.5;
        leaf_d[7] = 7.5;
        // leaf names
        vector<string> leaf_names;
        leaf_names.push_back("leaf 0");
        leaf_names.push_back("leaf 1");
        leaf_names.push_back("leaf 2");
        leaf_names.push_back("leaf 3");
        leaf_names.push_back("leaf 4");
        leaf_names.push_back("leaf 5");
        leaf_names.push_back("leaf 6");
        leaf_names.push_back("leaf 7");
        ntree_t ntree(nedge_set, leaf_d, leaf_names);
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
        test2();

        return 0;
}
