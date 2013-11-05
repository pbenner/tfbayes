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
        normal_fnode_t* n1 = new normal_fnode_t();
        normal_vnode_t* n2 = new normal_vnode_t();
         gamma_vnode_t* n3 = new  gamma_vnode_t();

        vector<variable_node_i*> vnodes;
        vector<  factor_node_i*> fnodes;

        n1->link(1, *n2);
        n1->link(2, *n3);

        fnodes.push_back(n1);
        vnodes.push_back(n2);
        vnodes.push_back(n3);

        factor_graph_t fg1(vnodes, fnodes);
        factor_graph_t fg2(fg1);
}
