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

#include <cassert>
#include <cstdlib>

#include <tfbayes/phylotree/phylotree-parsetree.hh>

using namespace std;

pt_parsetree_t::pt_parsetree_t(
        nodetype_t type,
        size_t n_children,
        void *data,
        ...)
        : data(data), type(type), n_children(n_children) {
        va_list az;

        children = (pt_parsetree_t**)calloc(n_children, sizeof(pt_parsetree_t *));

        va_start(az, data);
        for(size_t i = 0; i < n_children; i++) {
                children[i] = va_arg(az, pt_parsetree_t *);
        }
        va_end(az);
}

pt_parsetree_t::~pt_parsetree_t() {
        if (data) {
                free(data);
        }
        free(children);
}

void
pt_parsetree_t::destroy() {
        for(size_t i = 0; i < n_children; i++) {
                children[i]->destroy();
        }
        delete(this);
}

pt_node_t*
pt_parsetree_t::convert() const
{
        pt_node_t* node;

        switch (type)
        {
        case NODE_N:
                assert(children[2]->type == DISTANCE_N);
                node = new pt_node_t(*(double *)children[2]->data,
                                     children[0]->convert(),
                                     children[1]->convert());
                break;
        case LEAF_N:
                assert(children[0]->type == NAME_N);
                assert(children[1]->type == DISTANCE_N);
                node = new pt_leaf_t(-1,
                                    *(double *)children[1]->data,
                                     (char   *)children[0]->data);
                break;
        case ROOT_N:
                assert(children[0]->type == NODE_N ||
                       children[0]->type == LEAF_N);
                assert(children[1]->type == NODE_N ||
                       children[1]->type == LEAF_N);
                node = new pt_root_t(children[1]->convert(),
                                     children[0]->convert());
                break;
        default:
                assert(false);
                break;
        }

        return node;
}

#define TAB_WIDTH 4

static
void indent(ostream& o, size_t nesting)
{
        for(size_t i = 0; i < nesting*TAB_WIDTH; i++) {
                o << " ";
        }
}

static
ostream& print_parsetree(
        ostream& o,
        pt_parsetree_t* const tree,
        size_t nesting)
{
        indent(o, nesting);

        switch (tree->type)
        {
        case NAME_N:
                o << "NAME("
                  << (char *)tree->data
                  << ")"
                  << endl;
                break;
        case DISTANCE_N:
                o << "DISTANCE("
                  << *(double *)tree->data
                  << ")"
                  << endl;
                break;
        case NODE_N:
                o << "NODE"
                  << endl;
                break;
        case LEAF_N:
                o << "LEAF"
                  << endl;
                break;
        case ROOT_N:
                o << "ROOT"
                  << endl;
                break;
        default:
                break;
        }
        for (size_t i = 0; i < tree->n_children; i++) {
                print_parsetree(o, tree->children[i], nesting+1);
        }
        return o;
}

ostream& operator<< (ostream& o, pt_parsetree_t* const tree) {
        print_parsetree(o, tree, 0);
        return o;
}
