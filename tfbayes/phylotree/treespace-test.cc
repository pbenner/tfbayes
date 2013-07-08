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

#include <glpk.h>

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
        // print ntree
        cout << ntree
             << endl;
        // print exported tree
        pt_root_t* tree = ntree.export_tree();
        cout << "exporting tree:" << endl
             << tree              << endl;
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
        s4.insert(8);
        s4.insert(9);
        set<size_t> s5;
        s5.insert(0);
        s5.insert(1);
        s5.insert(8);
        s5.insert(9);
        set<size_t> s6;
        s6.insert(1);
        s6.insert(8);
        s6.insert(9);
        set<size_t> s7;
        s7.insert(8);
        s7.insert(9);
        nedge_t e1(9, s1, 1);
        nedge_t e2(9, s2, 1);
        nedge_t e3(9, s3, 1);
        nedge_t e4(9, s4, 1);
        nedge_t e5(9, s5, 1);
        nedge_t e6(9, s6, 1);
        nedge_t e7(9, s7, 1);
        // construct tree
        nedge_set_t nedge_set;
        // internal edge lengths
        nedge_set.push_back(e1);
        nedge_set.push_back(e2);
        nedge_set.push_back(e3);
        nedge_set.push_back(e4);
        nedge_set.push_back(e5);
        nedge_set.push_back(e6);
        nedge_set.push_back(e7);
        // leaf edge lengths
        vector<double> leaf_d(10, 1);
        leaf_d[0] = 0.5;
        leaf_d[1] = 1.5;
        leaf_d[2] = 2.5;
        leaf_d[3] = 3.5;
        leaf_d[4] = 4.5;
        leaf_d[5] = 5.5;
        leaf_d[6] = 6.5;
        leaf_d[7] = 7.5;
        leaf_d[8] = 8.5;
        leaf_d[9] = 9.5;
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
        leaf_names.push_back("leaf 8");
        leaf_names.push_back("leaf 9");
        ntree_t ntree(nedge_set, leaf_d, leaf_names);
        // print tree
        pt_root_t* tree = ntree.export_tree();
        cout << "exporting tree:" << endl
             << tree              << endl;
        tree->destroy();
}

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

void test5()
{
        cout << "test 5" << endl
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
        set<size_t> s4;
        s4.insert(2);
        s4.insert(3);
        set<size_t> s5;
        s5.insert(4);
        s5.insert(5);
        set<size_t> s6;
        s6.insert(0);
        s6.insert(1);
        nedge_t e1(5, s1, 5.5);
        nedge_t e2(5, s2, 2.5);
        nedge_t e3(5, s3, 1.5);
        nedge_t e4(5, s4, 1.5);
        nedge_t e5(5, s5, 1.5);
        nedge_t e6(5, s6, 1.5);
        // construct trees
        nedge_set_t nedge_set1;
        nedge_set_t nedge_set2;
        // internal edge lengths
        nedge_set1.push_back(e1);
        nedge_set1.push_back(e2);
        nedge_set1.push_back(e3);
        nedge_set2.push_back(e4);
        nedge_set2.push_back(e5);
        nedge_set2.push_back(e6);
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
        ntree_t nt1(nedge_set1, leaf_d, leaf_names);
        ntree_t nt2(nedge_set2, leaf_d, leaf_names);
        // print tree
        pt_root_t* t1 = nt1.export_tree();
        pt_root_t* t2 = nt2.export_tree();
        cout << "exporting tree:" << endl;
        cout << "exporting tree:" << endl
             << t1                << endl;
        cout << "exporting tree:" << endl
             << t2                << endl;
        t1->destroy();
        t2->destroy();

        geodesic_t geodesic(nt1, nt2);
        cout << "geodesic length: " << geodesic.length() << endl;
        pt_root_t* t3 = geodesic(1.0).export_tree();
        cout << "exporting tree:" << endl
             << t3                << endl;
        t3->destroy();
}

void test6()
{
        cout << "test 6" << endl
             <<"------------------------------------------------------------"
             << endl;
        // prepare splits
        set<size_t> s1;
        s1.insert(0);
        s1.insert(1);
        set<size_t> s2;
        s2.insert(2);
        s2.insert(3);
        set<size_t> s3;
        s3.insert(4);
        s3.insert(5);
        s3.insert(6);
        set<size_t> s4;
        s4.insert(0);
        s4.insert(2);
        s4.insert(3);
        set<size_t> s5;
        s5.insert(5);
        s5.insert(6);
        set<size_t> s6;
        s6.insert(4);
        s6.insert(5);
        nedge_t e1(6, s1, 5.5);
        nedge_t e2(6, s2, 2.5);
        nedge_t e31(6, s3, 1.5);
        nedge_t e32(6, s3, 2.5);
        nedge_t e4(6, s4, 0.5);
        nedge_t e5(6, s5, 3.5);
        nedge_t e6(6, s6, 4.5);
        // construct trees
        nedge_set_t nedge_set1;
        nedge_set_t nedge_set2;
        // internal edge lengths
        nedge_set1.push_back(e1);
        nedge_set1.push_back(e2);
        nedge_set1.push_back(e31);
        nedge_set1.push_back(e5);
        nedge_set2.push_back(e2);
        nedge_set2.push_back(e32);
        nedge_set2.push_back(e4);
        nedge_set2.push_back(e6);
        // leaf edge lengths
        vector<double> leaf_d(7, 1);
        leaf_d[0] = 0.5;
        leaf_d[1] = 1.5;
        leaf_d[2] = 2.5;
        leaf_d[3] = 3.5;
        leaf_d[4] = 4.5;
        leaf_d[5] = 5.5;
        leaf_d[6] = 6.5;
        // leaf names
        vector<string> leaf_names;
        leaf_names.push_back("leaf 0");
        leaf_names.push_back("leaf 1");
        leaf_names.push_back("leaf 2");
        leaf_names.push_back("leaf 3");
        leaf_names.push_back("leaf 4");
        leaf_names.push_back("leaf 5");
        leaf_names.push_back("leaf 6");
        ntree_t nt1(nedge_set1, leaf_d, leaf_names);
        ntree_t nt2(nedge_set2, leaf_d, leaf_names);
        // print tree
        pt_root_t* t1 = nt1.export_tree();
        pt_root_t* t2 = nt2.export_tree();
        cout << "exporting tree:" << endl
             << t1                << endl;
        cout << "exporting tree:" << endl
             << t2                << endl;
        t1->destroy();
        t2->destroy();

        cout << nt1 << endl;
        cout << nt2 << endl;

        geodesic_t geodesic(nt1, nt2);
        cout << "geodesic length: " << geodesic.length() << endl;
        pt_root_t* t3 = geodesic(1.0).export_tree();
        cout << "exporting tree:" << endl
             << t3                << endl;
        t3->destroy();
}

void test7()
{
        cout << "test 7" << endl
             <<"------------------------------------------------------------"
             << endl;
        // prepare splits
        set<size_t> s1;
        s1.insert(0);
        s1.insert(1);
        set<size_t> s2;
        s2.insert(1);
        s2.insert(2);
        set<size_t> s3;
        s3.insert(0);
        s3.insert(1);
        s3.insert(2);
        set<size_t> s4;
        s4.insert(4);
        s4.insert(5);
        set<size_t> s5;
        s5.insert(3);
        s5.insert(4);
        nedge_t e1(5, s1, 1.0);
        nedge_t e2(5, s2, 1.0);
        nedge_t e3(5, s3, 1.0);
        nedge_t e4(5, s4, 1.0);
        nedge_t e5(5, s5, 1.0);
        // construct trees
        nedge_set_t nedge_set1;
        nedge_set_t nedge_set2;
        // internal edge lengths
        nedge_set1.push_back(e1);
        nedge_set1.push_back(e3);
        nedge_set1.push_back(e4);
        nedge_set2.push_back(e2);
        nedge_set2.push_back(e3);
        nedge_set2.push_back(e5);
        // leaf edge lengths
        vector<double> leaf_d(6, 1);
        leaf_d[0] = 1.5;
        leaf_d[1] = 1.5;
        leaf_d[2] = 1.5;
        leaf_d[3] = 1.5;
        leaf_d[4] = 1.5;
        leaf_d[5] = 1.5;
        // leaf names
        vector<string> leaf_names;
        leaf_names.push_back("leaf 0");
        leaf_names.push_back("leaf 1");
        leaf_names.push_back("leaf 2");
        leaf_names.push_back("leaf 3");
        leaf_names.push_back("leaf 4");
        leaf_names.push_back("leaf 5");
        ntree_t nt1(nedge_set1, leaf_d, leaf_names);
        ntree_t nt2(nedge_set2, leaf_d, leaf_names);
        // print tree
        pt_root_t* t1 = nt1.export_tree();
        pt_root_t* t2 = nt2.export_tree();
        cout << "exporting tree:" << endl
             << t1                << endl;
        cout << "exporting tree:" << endl
             << t2                << endl;
        t1->destroy();
        t2->destroy();

        geodesic_t geodesic(nt1, nt2);
        cout << "npath list:" << endl
             << geodesic.npath_list();
        cout << "geodesic length: " << geodesic.length() << endl;
        pt_root_t* t3 = geodesic(1.0).export_tree();
        cout << "exporting tree:" << endl
             << t3                << endl;
        t3->destroy();
}

void test8()
{
        cout << "test 8" << endl
             << "------------------------------------------------------------"
             << endl;

        // Tree 1
        //////////////////////////////////////////////////////////////////////
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
        nedge_t e1(4, s1, 1.0, "e1");
        nedge_t e3_1(4, s3, 0.0, "e3");
        nedge_t e3_2(4, s3, 1.0, "e3");
        nedge_t e4(4, s4, 1.0, "e4");
        nedge_t e5(4, s5, 1.0, "e5");
        // construct trees
        nedge_set_t nedge_set1;
        nedge_set_t nedge_set2;
        // internal edge lengths
        nedge_set1.push_back(e1);
        //nedge_set1.push_back(e3_1);
        nedge_set2.push_back(e3_2);
        nedge_set2.push_back(e4);
        // leaf edge lengths
        vector<double> leaf_d(5, 1);
        leaf_d[0] = 1.5;
        leaf_d[1] = 1.5;
        leaf_d[2] = 1.5;
        leaf_d[3] = 1.5;
        leaf_d[4] = 1.5;
        // leaf names
        vector<string> leaf_names;
        leaf_names.push_back("leaf 0");
        leaf_names.push_back("leaf 1");
        leaf_names.push_back("leaf 2");
        leaf_names.push_back("leaf 3");
        leaf_names.push_back("leaf 4");
        ntree_t nt1(nedge_set1, leaf_d, leaf_names);
        ntree_t nt2(nedge_set2, leaf_d, leaf_names);

        geodesic_t geodesic(nt1, nt2);
        cout << geodesic.npath_list() << endl;
}

void test9()
{
        cout << "test 9" << endl
             << "------------------------------------------------------------"
             << endl;

        // Tree 1
        //////////////////////////////////////////////////////////////////////
        set<size_t> s1;
        s1.insert(1);
        s1.insert(2);
        set<size_t> s2;
        s2.insert(4);
        s2.insert(5);
        set<size_t> s3;
        s3.insert(0);
        s3.insert(4);
        s3.insert(5);
        set<size_t> s4;
        s4.insert(0);
        s4.insert(5);
        set<size_t> s5;
        s5.insert(0);
        s5.insert(3);
        nedge_t e1(5, s1, 1.0);
        nedge_t e2(5, s2, 1.0);
        nedge_t e3(5, s3, 1.0);
        nedge_t e4(5, s4, 1.0);
        nedge_t e5(5, s5, 0.0);
        // construct trees
        nedge_set_t nedge_set1;
        nedge_set_t nedge_set2;
        // internal edge lengths
        nedge_set1.push_back(e1);
        nedge_set1.push_back(e2);
        //nedge_set1.push_back(e5);
        nedge_set2.push_back(e1);
        nedge_set2.push_back(e3);
        nedge_set2.push_back(e4);
        // leaf edge lengths
        vector<double> leaf_d(6, 1);
        leaf_d[0] = 1.5;
        leaf_d[1] = 1.5;
        leaf_d[2] = 1.5;
        leaf_d[3] = 1.5;
        leaf_d[4] = 1.5;
        leaf_d[5] = 1.5;
        // leaf names
        vector<string> leaf_names;
        leaf_names.push_back("leaf 0");
        leaf_names.push_back("leaf 1");
        leaf_names.push_back("leaf 2");
        leaf_names.push_back("leaf 3");
        leaf_names.push_back("leaf 4");
        leaf_names.push_back("leaf 5");
        ntree_t nt1(nedge_set1, leaf_d, leaf_names);
        ntree_t nt2(nedge_set2, leaf_d, leaf_names);

        geodesic_t geodesic(nt1, nt2);
        cout << geodesic.npath_list() << endl;

        // print length
        cout << "length: "
             << geodesic.length()
             << endl;

}

void test10()
{
        cout << "test 10" << endl
             << "------------------------------------------------------------"
             << endl;

        // Tree 1
        //////////////////////////////////////////////////////////////////////
        set<size_t> s1;
        s1.insert(4);
        s1.insert(5);
        set<size_t> s2;
        s2.insert(2);
        s2.insert(3);
        set<size_t> s3;
        s3.insert(1);
        s3.insert(2);
        s3.insert(3);
        nedge_t e1(5, s1, 1.0);
        nedge_t e2(5, s2, 1.0);
        nedge_t e3(5, s3, 1.0);
        // construct trees
        nedge_set_t nedge_set1;
        nedge_set_t nedge_set2;
        // internal edge lengths
        nedge_set1.push_back(e1);
        nedge_set1.push_back(e2);
        nedge_set1.push_back(e3);
        //nedge_set1.push_back(e5);
        nedge_set2.push_back(e2);
        nedge_set2.push_back(e1);
        nedge_set2.push_back(e3);
        // leaf edge lengths
        vector<double> leaf_d(6, 1);
        leaf_d[0] = 1.5;
        leaf_d[1] = 1.5;
        leaf_d[2] = 1.5;
        leaf_d[3] = 1.5;
        leaf_d[4] = 1.5;
        leaf_d[5] = 1.5;
        // leaf names
        vector<string> leaf_names;
        leaf_names.push_back("mm9");
        leaf_names.push_back("papHam1");
        leaf_names.push_back("gorGor1");
        leaf_names.push_back("panTro2");
        leaf_names.push_back("ponAbe2");
        leaf_names.push_back("hg19");
        ntree_t nt1(nedge_set1, leaf_d, leaf_names);
        ntree_t nt2(nedge_set2, leaf_d, leaf_names);

        pt_root_t* t1 = nt1.export_tree();
        cout << "exporting tree:" << endl
             << t1                << endl;
        t1->destroy();

        pt_root_t* t2 = nt2.export_tree();
        cout << "exporting tree:" << endl
             << t2                << endl;
        t2->destroy();
}

void test11()
{
        cout << "test 11" << endl
             << "------------------------------------------------------------"
             << endl;

        // Tree 1
        //////////////////////////////////////////////////////////////////////
        set<size_t> s1;
        s1.insert(1);
        s1.insert(2);
        set<size_t> s2;
        s2.insert(1);
        s2.insert(2);
        s2.insert(3);
        set<size_t> s3;
        s3.insert(1);
        s3.insert(2);
        s3.insert(3);
        s3.insert(4);
        set<size_t> s4;
        s4.insert(1);
        s4.insert(2);
        s4.insert(3);
        s4.insert(4);
        s4.insert(5);
        set<size_t> s5;
        s5.insert(1);
        s5.insert(2);
        s5.insert(3);
        s5.insert(4);
        s5.insert(5);
        s5.insert(6);
        set<size_t> s6;
        s6.insert(7);
        s6.insert(8);
        set<size_t> s7;
        s7.insert(7);
        s7.insert(8);
        s7.insert(9);
        set<size_t> s8;
        s8.insert(10);
        s8.insert(11);
        set<size_t> s9;
        s9.insert(10);
        s9.insert(11);
        s9.insert(12);
        set<size_t> s10;
        s10.insert(7);
        s10.insert(8);
        s10.insert(9);
        s10.insert(10);
        s10.insert(11);
        s10.insert(12);
        nedge_t e1 (12,  s1, 1.0);
        nedge_t e2 (12,  s2, 1.0);
        nedge_t e3 (12,  s3, 1.0);
        nedge_t e4 (12,  s4, 1.0);
        nedge_t e5 (12,  s5, 1.0);
        nedge_t e6 (12,  s6, 1.0);
        nedge_t e7 (12,  s7, 1.0);
        nedge_t e8 (12,  s8, 1.0);
        nedge_t e9 (12,  s9, 1.0);
        nedge_t e10(12, s10, 1.0);
        // construct tree
        nedge_set_t nedge_set;
        // internal edge lengths
        nedge_set.push_back(e1);
        nedge_set.push_back(e2);
        nedge_set.push_back(e3);
        nedge_set.push_back(e4);
        nedge_set.push_back(e5);
        nedge_set.push_back(e6);
        nedge_set.push_back(e7);
        nedge_set.push_back(e8);
        nedge_set.push_back(e9);
        nedge_set.push_back(e10);
        // leaf edge lengths
        vector<double> leaf_d(13, 1);
        // leaf names
        vector<string> leaf_names;
        leaf_names.push_back("ochPri2");
        leaf_names.push_back("hg19");
        leaf_names.push_back("panTro2");
        leaf_names.push_back("gorGor1");
        leaf_names.push_back("ponAbe2");
        leaf_names.push_back("papHam1");
        leaf_names.push_back("calJac1");
        leaf_names.push_back("speTri1");
        leaf_names.push_back("oryCun2");
        leaf_names.push_back("cavPor3");
        leaf_names.push_back("mm9");
        leaf_names.push_back("rn4");
        leaf_names.push_back("dipOrd1");
        // tree
        ntree_t nt(nedge_set, leaf_d, leaf_names);
        cout << "tree:" << endl
             << nt      << endl;

        pt_root_t* tree = nt.export_tree();
        cout << "exporting tree:" << endl
             << tree              << endl;
        tree->destroy();
}

void test12()
{
        cout << "test 12" << endl
             << "------------------------------------------------------------"
             << endl;

        //////////////////////////////////////////////////////////////////////
        set<size_t> s1;
        s1.insert(1);
        s1.insert(2);
        set<size_t> s2;
        s2.insert(4);
        s2.insert(5);
        set<size_t> s3;
        s3.insert(6);
        s3.insert(7);
        set<size_t> s4;
        s4.insert(4);
        s4.insert(5);
        s4.insert(6);
        s4.insert(7);
        nedge_t e1 (7,  s1, 1.0);
        nedge_t e2 (7,  s2, 1.0);
        nedge_t e3 (7,  s3, 1.0);
        nedge_t e4 (7,  s4, 1.0);
        // construct tree
        nedge_set_t nedge_set;
        // internal edge lengths
        nedge_set.push_back(e1);
        nedge_set.push_back(e2);
        nedge_set.push_back(e3);
        nedge_set.push_back(e4);
        // leaf edge lengths
        vector<double> leaf_d(8, 1);
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
        // tree
        ntree_t nt(nedge_set, leaf_d, leaf_names);
        cout << "tree:" << endl
             << nt      << endl;

        pt_root_t* tree = nt.export_tree();
        cout << "exporting tree:" << endl
             << tree              << endl;
        tree->destroy();
}

int
main(void)
{
        test1(); cout << endl;
        test2(); cout << endl;
        test3(); cout << endl;
        test4(); cout << endl;
        test5(); cout << endl;
        test6(); cout << endl;
        test7(); cout << endl;
        test8(); cout << endl;
        test9(); cout << endl;
        test10(); cout << endl;
        test11(); cout << endl;
        test12(); cout << endl;
        glp_free_env();

        return 0;
}
