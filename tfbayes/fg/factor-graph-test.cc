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
        factor_set_t fnodes;
        variable_set_t vnodes;

        vector<double> data(1, 1.42);

        fnodes += new normal_fnode_t("f1", 1, 2);
        vnodes += new normal_data_t ("v1");
        fnodes += new normal_fnode_t("f2", 2, 2);
        vnodes += new normal_vnode_t("v2");
        fnodes += new gamma_fnode_t ("f3", 1, 2);
        vnodes += new gamma_vnode_t ("v3");

        fnodes["f1"]->link("output",    *vnodes["v1"]);
        fnodes["f1"]->link("mean",      *vnodes["v2"]);
        fnodes["f1"]->link("precision", *vnodes["v3"]);
        fnodes["f2"]->link("output",    *vnodes["v2"]);
        fnodes["f3"]->link("output",    *vnodes["v3"]);

        vnodes["v1"]->condition(data);

        factor_graph_t fg1(fnodes, vnodes, 1);
        factor_graph_t fg2(fg1);

        fg1();
}

void
test2()
{
        factor_graph_t fg;

        vector<double> data;

        data.push_back(1.42);
        data.push_back(1.81);

        fg += normal_fnode_t("f1", 1, 2, 2);
        fg += normal_data_t ("v1");
        fg += normal_fnode_t("f2", 2, 2);
        fg += normal_vnode_t("v2");
        fg += gamma_fnode_t ("f3", 1, 2);
        fg += gamma_vnode_t ("v3");

        fg.link("f1", "output",    "v1");
        fg.link("f1", "mean",      "v2");
        fg.link("f1", "precision", "v3");
        fg.link("f2", "output",    "v2");
        fg.link("f3", "output",    "v3");

        fg.variable_node("v1")->condition(data);

        factor_graph_t fg2(fg);

        fg();

        cout << "mean: " << fg.distribution("v2")->moments()[0] << endl
             << "mean: " << fg.distribution("v3")->moments()[0] << endl;
}

void
test3()
{
        factor_graph_t fg;

        vector<double> data1(1, 1.42);
        vector<double> data2(1, 1.81);

        fg += normal_fnode_t("f1", 1, 2);
        fg += normal_data_t ("v1");
        fg.link("f1", "output", "v1");
        fg.replicate(1);

        fg.variable_node("v1", 0)->condition(data1);
        fg.variable_node("v1", 1)->condition(data2);

        fg += normal_fnode_t("f2", 2, 2);
        fg += normal_vnode_t("v2");
        fg += gamma_fnode_t ("f3", 1, 2);
        fg += gamma_vnode_t ("v3");

        fg.link("f1", "mean",      "v2");
        fg.link("f1", "precision", "v3");
        fg.link("f2", "output",    "v2");
        fg.link("f3", "output",    "v3");

        factor_graph_t fg2(fg);

        fg();

        cout << "mean: " << fg.distribution("v2")->moments()[0] << endl
             << "mean: " << fg.distribution("v3")->moments()[0] << endl;
}

int
main()
{
        test3();
}
