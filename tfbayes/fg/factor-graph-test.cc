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

#include <iostream>
#include <vector>

#include <factor-graph.hh>
#include <variational.hh>

using namespace std;

void
test1()
{
        factor_graph_t fg1;

        matrix<double> data(2, 1);
        data[0][0] = 1.42;
        data[1][0] = 1.81;

        fg1 += normal_fnode_t("f1", 1, 2);
        fg1 += normal_vnode_t("v1");
        fg1 += normal_fnode_t("f2", 2, 2);
        fg1 += normal_vnode_t("v2");
        fg1 += gamma_fnode_t ("f3", 1, 2);
        fg1 += gamma_vnode_t ("v3");

        fg1.link("f1:output",    "v1");
        fg1.link("f1:mean",      "v2");
        fg1.link("f1:precision", "v3");
        fg1.link("f2:output",    "v2");
        fg1.link("f3:output",    "v3");

        fg1.variable_node("v1")->condition(data);

        factor_graph_t fg2 = fg1;

        fg2(5);

        cout << "mean: " << fg2.distribution("v2")->moments()[0] << endl
             << "mean: " << fg2.distribution("v3")->moments()[0] << endl;
}

void
test2()
{
        factor_graph_t fg;

        matrix<double> data1(1, 1);
        data1[0][0] = 1.42;
        matrix<double> data2(1, 1);
        data2[0][0] = 1.81;

        fg += normal_fnode_t("f1", 1, 2);
        fg += normal_vnode_t("v1");
        fg.link("f1:output", "v1");
        fg.replicate(1);

        fg.variable_node("v1", 0)->condition(data1);
        fg.variable_node("v1", 1)->condition(data2);

        fg += normal_fnode_t("f2", 2, 2);
        fg += normal_vnode_t("v2");
        fg += gamma_fnode_t ("f3", 1, 2);
        fg += gamma_vnode_t ("v3");

        fg.link("f1:mean",      "v2");
        fg.link("f1:precision", "v3");
        fg.link("f2:output",    "v2");
        fg.link("f3:output",    "v3");

        factor_graph_t fg2(fg);

        fg();

        cout << "mean: " << fg.distribution("v2")->moments()[0] << endl
             << "mean: " << fg.distribution("v3")->moments()[0] << endl;
}

void
test3()
{
        factor_graph_t fg;

        matrix<double> data(5,1);
        data[0][0] = 1;
        data[1][0] = 1;
        data[2][0] = 2;
        data[3][0] = 0;
        data[4][0] = 1;

        vector<double> theta(3, 0.0);
        vector<double> alpha(3, 1.0);

        fg += categorical_fnode_t("f1", theta);
        fg += categorical_vnode_t("v1", 3);
        fg += dirichlet_fnode_t("f2", alpha);
        fg += dirichlet_vnode_t("v2", 3);

        fg.link("f1:output",  "v1");
        fg.link("f1:theta",   "v2");
        fg.link("f2:output",  "v2");

        fg.variable_node("v1", 0)->condition(data);

        fg();
}

void
test4()
{
        factor_graph_t fg1;

        matrix<double> data(1,1);
        data[0][0] = 0.0;

        // pseudocounts for the dirichlet prior
        vector<double> alpha(2, 1.0);

        mixture_fnode_t m("mixture_fnode_1");
        m += normal_fnode_t("normal_1");
        m += normal_fnode_t("normal_2");

        fg1 += m;
        fg1 += normal_vnode_t("x");
        fg1 += normal_vnode_t("normal_vnode_1");
        fg1 += normal_vnode_t("normal_vnode_2");
        fg1 += normal_fnode_t("normal_fnode_1", 0, 10);
        fg1 += normal_fnode_t("normal_fnode_2", 0, 10);
        fg1 += categorical_fnode_t("categorical_fnode_1", 2);
        fg1 += categorical_vnode_t("categorical_vnode_1", 2);
        fg1 += dirichlet_fnode_t("dirichlet_fnode_1", alpha);
        fg1 += dirichlet_vnode_t("dirichlet_vnode_1", 2);

        fg1.link("mixture_fnode_1:output",        "x");
        fg1.link("mixture_fnode_1:normal_1:mean", "normal_vnode_1");
        fg1.link("mixture_fnode_1:normal_2:mean", "normal_vnode_2");
        fg1.link("mixture_fnode_1:indicator",     "categorical_vnode_1");
        fg1.link("categorical_fnode_1:output",    "categorical_vnode_1");
        fg1.link("categorical_fnode_1:theta",     "dirichlet_vnode_1");
        fg1.link("dirichlet_fnode_1:output",      "dirichlet_vnode_1");
        fg1.link("normal_fnode_1:output",         "normal_vnode_1");
        fg1.link("normal_fnode_2:output",         "normal_vnode_2");

        fg1.variable_node("x", 0)->condition(data);

        factor_graph_t fg2 = fg1;

        vector<double> bound = fg2();

        for (size_t i = 0; i < bound.size(); i++) {
                cout << "bound: " << bound[i] << endl;
        }

        cout << "theta: " << fg2.distribution("categorical_vnode_1")->moments()[0] << endl
             << "theta: " << fg2.distribution("categorical_vnode_1")->moments()[1] << endl;
        cout << "alpha: " << fg2.distribution("dirichlet_vnode_1")->parameters()[0]+1.0 << endl
             << "alpha: " << fg2.distribution("dirichlet_vnode_1")->parameters()[1]+1.0 << endl;
}

int
main()
{
        test4();
}
