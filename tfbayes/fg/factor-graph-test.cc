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

int
main()
{
        factor_set_t fnodes;
        variable_set_t vnodes;

        fnodes += new normal_fnode_t(1,2, "f1");
        vnodes += new data_vnode_t(1.42, "v1");

        fnodes += new normal_fnode_t(2,2, "f2");
        vnodes += new normal_vnode_t("v2");

        fnodes += new gamma_fnode_t(1,2, "f3");
        vnodes += new gamma_vnode_t("v3");

        fnodes["f1"]->link("output",    *vnodes["v1"]);
        fnodes["f1"]->link("mean",      *vnodes["v2"]);
        fnodes["f1"]->link("precision", *vnodes["v3"]);
        fnodes["f2"]->link("output",    *vnodes["v2"]);
        fnodes["f3"]->link("output",    *vnodes["v3"]);

        std::cout << std::endl;
        factor_graph_t fg1(fnodes, vnodes, 1);
        //factor_graph_t fg2(fg1);

        fg1();

        // cout << "mean: " << (*v2)().moment<1>() << endl
        //      << "mean: " << (*v3)().moment<1>() << endl;
}
