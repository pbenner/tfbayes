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
        cout << "test 1" << endl
             <<"------------------------------------------------------------"
             << endl;
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
        cout << "test 2" << endl
             <<"------------------------------------------------------------"
             << endl;
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

#include <glpk.h>

void test3()
{
        cout << "test 3" << endl
             <<"------------------------------------------------------------"
             << endl;
        glp_prob *lp;
        glp_smcp parm;
        int ia[1+1000], ja[1+1000];
        double ar[1+1000], z, x1, x2, x3, x4;
        lp = glp_create_prob();
        glp_init_smcp(&parm);
        glp_set_obj_dir(lp, GLP_MIN);
        //
        parm.msg_lev = GLP_MSG_OFF;
        //
        glp_add_rows(lp, 3);
        glp_add_cols(lp, 4);
        //
        glp_set_row_bnds(lp, 1, GLP_LO, 1.0, 0.0);
        glp_set_row_bnds(lp, 2, GLP_LO, 1.0, 0.0);
        glp_set_row_bnds(lp, 3, GLP_LO, 1.0, 0.0);
        //
        glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 0.0);
        glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0);
        glp_set_col_bnds(lp, 3, GLP_LO, 0.0, 0.0);
        glp_set_col_bnds(lp, 4, GLP_LO, 0.0, 0.0);
        //
        glp_set_col_kind(lp, 1, GLP_IV);
        glp_set_col_kind(lp, 2, GLP_IV);
        glp_set_col_kind(lp, 3, GLP_IV);
        glp_set_col_kind(lp, 4, GLP_IV);
        //
        glp_set_obj_coef(lp, 1, 2.0);
        glp_set_obj_coef(lp, 2, 1.0);
        glp_set_obj_coef(lp, 3, 1.0);
        glp_set_obj_coef(lp, 4, 1.0);
        //
        ia[ 1] = 1, ja[ 1] = 1, ar[ 1] = 1.0;
        ia[ 2] = 1, ja[ 2] = 3, ar[ 2] = 1.0;
        ia[ 3] = 2, ja[ 3] = 2, ar[ 3] = 1.0;
        ia[ 4] = 2, ja[ 4] = 3, ar[ 4] = 1.0;
        ia[ 5] = 3, ja[ 5] = 2, ar[ 5] = 1.0;
        ia[ 6] = 3, ja[ 6] = 4, ar[ 6] = 1.0;
        glp_load_matrix(lp, 6, ia, ja, ar);
        //
        glp_simplex(lp, &parm);
        //
        if (glp_get_status(lp) == GLP_OPT) {
                printf("Optimal solution found.\n");
        }
        //
        z = glp_get_obj_val(lp);
        x1 = glp_get_col_prim(lp, 1);
        x2 = glp_get_col_prim(lp, 2);
        x3 = glp_get_col_prim(lp, 3);
        x4 = glp_get_col_prim(lp, 4);
        printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g; x4 = %g\n",
               z, x1, x2, x3, x4);
        glp_delete_prob(lp);
}

void test4()
{
        cout << "test 4" << endl
             <<"------------------------------------------------------------"
             << endl;
        // prepare splits
        set<size_t> s1;
        s1.insert(1);
        s1.insert(2);
        set<size_t> s3;
        s3.insert(0);
        s3.insert(4);
        set<size_t> s4;
        s4.insert(2);
        s4.insert(3);
        set<size_t> s5;
        s5.insert(0);
        s5.insert(1);
        nedge_t e1(4, s1, 0.5);
        nedge_t e3(4, s3, 1.0);
        nedge_t e4(4, s4, 1.0);
        nedge_t e5(4, s5, 0.5);
        // construct edge sets
        nedge_set_t nedge_set1;
        nedge_set_t nedge_set2;
        nedge_set1.push_back(e1);
        nedge_set1.push_back(e3);
        nedge_set2.push_back(e4);
        nedge_set2.push_back(e5);

        incompatibility_graph_t igraph(nedge_set1, nedge_set2);
        vertex_cover_t vertex_cover = igraph.min_weight_cover();

        cout << "weight: " << vertex_cover.weight << endl;
        cout << "set a:" << endl;
        for (size_t i = 0; i < vertex_cover.a.size(); i++) {
                cout << vertex_cover.a[i]
                     << endl;
        }
        cout << "set a_comp:" << endl;
        for (size_t i = 0; i < vertex_cover.a_comp.size(); i++) {
                cout << vertex_cover.a_comp[i]
                     << endl;
        }
        cout << "set b:" << endl;
        for (size_t i = 0; i < vertex_cover.b.size(); i++) {
                cout << vertex_cover.b[i]
                     << endl;
        }
        cout << "set b_comp:" << endl;
        for (size_t i = 0; i < vertex_cover.b_comp.size(); i++) {
                cout << vertex_cover.b_comp[i]
                     << endl;
        }
}

int
main(void)
{
        test1(); cout << endl;
        test2(); cout << endl;
        test3(); cout << endl;
        test4(); cout << endl;

        return 0;
}
