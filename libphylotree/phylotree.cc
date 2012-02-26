/* Copyright (C) 2012 Philipp Benner
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

#include <phylotree.hh>

using namespace std;

#define TAB_WIDTH 2

static
void indent(ostream& o, size_t nesting)
{
        for(size_t i = 0; i < nesting*TAB_WIDTH; i++) {
                o << " ";
        }
}

static
ostream& print_phylotree(
        ostream& o,
        pt_node_t* const node,
        size_t nesting)
{
        if (node->root()) {
                cout << "(root";
                print_phylotree(o, node->left,  nesting+1);
                print_phylotree(o, node->right, nesting+1);
                cout << ")"
                     << endl;
        }
        else if (node->leaf()) {
                cout << endl;
                indent(o, nesting);
                cout << "("
                     << node->name
                     << " "
                     << node->d;
                if (node->x != -1) {
                        cout << " "
                             << node->x;
                }
                cout << ")";
        }
        else {
                cout << endl;
                indent(o, nesting);
                cout << "(node "
                     << node->d;
                print_phylotree(o, node->left,  nesting+1);
                print_phylotree(o, node->right, nesting+1);
                cout << ")";
        }

        return o;
}

ostream& operator<< (ostream& o, pt_node_t* const node) {
        print_phylotree(o, node, 0);
        return o;
}
