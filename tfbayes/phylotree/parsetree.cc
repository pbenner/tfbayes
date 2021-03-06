/* Copyright (C) 2012-2013 Philipp Benner
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
#include <new>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <boost/format.hpp>

#include <tfbayes/phylotree/parser.hh>
#include <tfbayes/phylotree/parsetree.hh>

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
        // non-recursive destruction of the tree list
        if (type != TREE_LIST_N) {
                for(size_t i = 0; i < n_children; i++) {
                        delete(children[i]);
                }
        }
        if (data) {
                free(data);
        }
        free(children);
}

void delete_tree_list(pt_parsetree_t* tree) {
        assert(tree->type == TREE_LIST_N);

        pt_parsetree_t *tmp;

        while (tree) {
                if (tree->n_children == 2) {
                        // delete phylogenetic tree
                        delete(tree->children[1]);
                        // go down in the AST
                        tmp  = tree;
                        tree = tree->children[0];
                        // free current node
                        delete(tmp);
                }
                else {
                        // delete phylogenetic tree
                        delete(tree->children[0]);
                        // free current node
                        delete(tree);
                        // end recursion
                        tree = NULL;
                }
        }
}

list<pt_root_t>
pt_parsetree_t::convert(size_t drop, size_t skip,
                        boost::optional<const pt_root_t&> ref_tree) const
{
        size_t n = 0;
        list<pt_root_t> tree_list;
        const pt_parsetree_t* pt = this;

        // find the number of trees in the list
        for (size_t i = 0; pt != NULL; i++) {
                if (pt->n_children == 2) {
                        pt = pt->children[0]; n++;
                }
                else {
                        pt = NULL;
                }
        }
        pt = this;

        // iterate throught the tree list non-recursively
        // since it might get very big!
        for (size_t i = 0; pt != NULL; i++) {
                assert(pt->type == TREE_LIST_N);

                if (pt->n_children == 2) {
                        assert(pt->children[1]->type == TREE_N);
                        // convert the tree on the right which is a
                        // phylogenetic tree
                        if (n-i >= drop && i % skip == 0) {
                                tree_list.push_front(pt->children[1]->convert(tree_list, ref_tree));
                        }
                        // on the left is the tree list
                        assert(pt->children[0]->type == TREE_LIST_N);
                        pt = pt->children[0];
                }
                else {
                        assert(pt->children[0]->type == TREE_N);
                        // we only have the phylogenetic tree
                        if (n-i >= drop && i % skip == 0) {
                                tree_list.push_front(pt->children[0]->convert(tree_list, ref_tree));
                        }
                        // end loop
                        pt = NULL;
                }
        }
        return tree_list;
}

pt_root_t
pt_parsetree_t::convert(const list<pt_root_t>& tree_list,
                        boost::optional<const pt_root_t&> ref_tree) const
{
        assert(type == TREE_N);
        assert(n_children == 1);
        assert(children[0]->type == NODE_LIST_N);
        assert(children[0]->n_children == 2);
        // make sure that there are at least three childs
        // attached to the root (including outgroup)
        pt_leaf_t* outgroup;
        pt_node_t* node;
        // if the right branch is a leaf use it as the outgroup
        if (children[0]->children[1]->type == LEAF_N) {
                outgroup = static_cast<pt_leaf_t*>(
                        children[0]->children[1]->convert());
                node = children[0]->children[0]->convert();
        }
        else {
                outgroup = NULL;
                node = children[0]->convert();
        }
        if (ref_tree) {
                return pt_root_t(*node, outgroup, *ref_tree);
        }
        else if (tree_list.begin() != tree_list.end()) {
                return pt_root_t(*node, outgroup, *tree_list.begin());
        }
        else {
                return pt_root_t(*node, outgroup);
        }
}

pt_node_t*
pt_parsetree_t::convert() const
{
        pt_node_t* child_left;
        pt_node_t* child_right;
        pt_node_t* node = NULL;

        switch (type)
        {
        case NODE_LIST_N:
                if (n_children == 1) {
                        return children[0]->convert();
                }
                else {
                        child_left  = children[0]->convert();
                        child_right = children[1]->convert();
                        try {
                                node = new pt_node_t(0.0, child_left, child_right);
                        }
                        catch (bad_alloc&) {
                                cerr << "Memory allocation failed!"
                                     << endl;
                        }
                }
                break;
        case NODE_N:
                assert(n_children == 2);
                assert(children[0]->type == NODE_LIST_N);
                assert(children[1]->type == DISTANCE_N);
                node    = children[0]->convert();
                node->d = *(double *)children[1]->data;
                break;
        case LEAF_N:
                assert(children[0]->type == NAME_N);
                assert(children[1]->type == DISTANCE_N);
                try {
                        node = new pt_leaf_t(*(double *)children[1]->data,
                                              (char   *)children[0]->data);
                }
                catch (bad_alloc&) {
                        cerr << "Memory allocation failed!"
                             << endl;
                }
                break;
        case TREE_N:
        case TREE_LIST_N:
                // this shouldn't happen
                assert(type != TREE_N);
                assert(type != TREE_LIST_N);
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
        case NODE_LIST_N:
                o << "NODE_LIST"
                  << endl;
                break;
        case LEAF_N:
                o << "LEAF"
                  << endl;
                break;
        case TREE_N:
                o << "TREE"
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

#include <cerrno>
#include <cstring>

int yylex_init   (yyscan_t* scanner);
int yylex_destroy(yyscan_t  scanner);
void yylex_set_input(yyscan_t scanner, FILE* file);

list<pt_root_t> parse_tree_list(FILE * file, size_t drop, size_t skip,
                                boost::optional<const pt_root_t&> ref_tree)
{
        list<pt_root_t> tree_list;
        context_t context;

        // initialize lexer
        yylex_init(&context.scanner);

        // set input file if necessary
        if (file != NULL) {
                yylex_set_input(context.scanner, file);
        }

        // parse input
        if (yyparse(&context)) {
                // error parsing
                return tree_list;
        }
        // convert AST to phylotree
        tree_list = context.pt_parsetree->convert(drop, skip, ref_tree);
 
        // free lexer memory
        yylex_destroy(context.scanner);

        // free AST
        delete_tree_list(context.pt_parsetree);
 
        return tree_list;
}

list<pt_root_t> parse_tree_list(const string& filename, size_t drop, size_t skip,
                                boost::optional<const pt_root_t&> ref_tree)
{
        FILE* yyin = NULL;
        if (filename != "") {
                yyin = fopen(filename.c_str(), "r");
                if (yyin == NULL) {
                        cerr << boost::format("Could not open tree file `%s': %s") % filename % strerror(errno)
                             << endl;
                        exit(EXIT_FAILURE);
                }
        }
        list<pt_root_t> tree_list = parse_tree_list(yyin, drop, skip, ref_tree);
        if (yyin) fclose(yyin);

        return tree_list;
}

// ostream
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, pt_parsetree_t* const tree) {
        print_parsetree(o, tree, 0);
        return o;
}
