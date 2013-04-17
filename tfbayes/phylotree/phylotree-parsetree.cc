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
#include <cassert>
#include <cstdlib>

#include <tfbayes/phylotree/phylotree-parser.h>
#include <tfbayes/phylotree/phylotree-parser.hh>
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

list<pt_root_t*>
pt_parsetree_t::convert() const
{
        list<pt_root_t*> tree_list;

        convert(tree_list);

        return tree_list;
}

pt_node_t*
pt_parsetree_t::convert(list<pt_root_t*>& tree_list) const
{
        pt_node_t* node = NULL;
        pt_root_t* root;

        switch (type)
        {
        case NODE_N:
                assert(children[2]->type == DISTANCE_N);
                node = new pt_node_t(*(double *)children[2]->data,
                                     children[0]->convert(tree_list),
                                     children[1]->convert(tree_list));
                break;
        case LEAF_N:
                assert(children[0]->type == NAME_N);
                assert(children[1]->type == DISTANCE_N);
                node = new pt_leaf_t(*(double *)children[1]->data,
                                      (char   *)children[0]->data);
                break;
        case ROOT_N:
                assert(children[0]->type == NODE_N ||
                       children[0]->type == LEAF_N);
                assert(children[1]->type == NODE_N ||
                       children[1]->type == LEAF_N);
                assert((n_children == 3 && children[2]->type == LEAF_N) ||
                        n_children == 2);
                node = new pt_root_t(children[0]->convert(tree_list),
                                     children[1]->convert(tree_list));
                if (n_children == 3) {
                        node->d = *(double *)children[2]->children[1]->data;
                        static_cast<pt_root_t*>(node)->outgroup_name =
                                string((char   *)children[2]->children[0]->data);
                }
                break;
        case TREE_N:
                assert(children[0]->type == ROOT_N);
                root = static_cast<pt_root_t*>(children[0]->convert(tree_list));
                tree_list.push_back(root);
                if (n_children == 2) {
                        assert(children[1]->type == TREE_N);
                        children[1]->convert(tree_list);
                }
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

// bison/flex interface
////////////////////////////////////////////////////////////////////////////////

typedef void* yyscan_t;
int yylex_init   (yyscan_t* scanner);
int yylex_destroy(yyscan_t  scanner);
void yylex_set_input(yyscan_t scanner, FILE* file);

list<pt_root_t*> parse_tree_list(FILE * file)
{
        list<pt_root_t*> tree_list;
        yyscan_t scanner;

        // initialize lexer
        yylex_init(&scanner);

        // set input file if necessary
        if (file != NULL) {
                yylex_set_input(scanner, file);
        }

        // parse input
        if (yyparse(scanner)) {
                // error parsing
                return tree_list;
        }
        // convert AST to phylotree
        tree_list = pt_parsetree->convert();
 
        // free lexer memory
        yylex_destroy(scanner);

        // free AST
        pt_parsetree->destroy();
 
        return tree_list;
}

// ostream
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, pt_parsetree_t* const tree) {
        print_parsetree(o, tree, 0);
        return o;
}
