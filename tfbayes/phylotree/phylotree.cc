/* Copyright (C) 2012, 2013 Philipp Benner
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

// pt_node_t
////////////////////////////////////////////////////////////////////////////////

pt_node_t::pt_node_t(double d,
                     pt_node_t* left,
                     pt_node_t* right,
                     const string& name,
                     bool is_root,
                     id_t id)
        : d(d), left(left), right(right), name(name),
          is_root(is_root), id(id)
{
        assert ((left == NULL && right == NULL) ||
                (left != NULL && right != NULL));
        
        if (leaf()) {
                n_leafs = 1;
                n_nodes = 1;
        }
        else {
                n_leafs = left->n_leafs + right->n_leafs;
                n_nodes = left->n_nodes + right->n_nodes + 1;
        }
}

pt_node_t::pt_node_t(const pt_node_t& node)
        : d(node.d), left(NULL), right(NULL),
          name(node.name), is_root(node.is_root),
          n_leafs(node.n_leafs)
{
        
        if (!node.leaf()) {
                left  = node.left ->clone();
                right = node.right->clone();
        }
}

pt_node_t*
pt_node_t::clone() const
{
        return new pt_node_t(*this);
}

void
pt_node_t::destroy()
{
        if (!leaf()) {
                left ->destroy();
                right->destroy();
        }
        delete(this);
}

bool
pt_node_t::leaf() const
{
        return left == NULL && right == NULL;
}

bool
pt_node_t::root() const
{
        return is_root;
}

double
pt_node_t::mutation_probability() const
{
        return 1.0-exp(-d);
}

void
pt_node_t::scale(double c)
{
        d *= c;
        if (!leaf()) {
                left->scale(c);
                right->scale(c);
        }
}

void
pt_node_t::get_leafs(leafs_t& leafs)
{
        if (!leaf()) {
                left ->get_leafs(leafs);
                right->get_leafs(leafs);
        }
}

void
pt_node_t::get_nodes(nodes_t& nodes)
{
        if (!root()) {
                nodes.insert(this);
        }
        if (!leaf()) {
                left ->get_nodes(nodes);
                right->get_nodes(nodes);
        }
}

// pt_leaf_t
////////////////////////////////////////////////////////////////////////////////

pt_leaf_t::pt_leaf_t(short x, double d, const std::string name)
        : pt_node_t(d, NULL, NULL, name),
          x(x)
{ }

void
pt_leaf_t::get_leafs(leafs_t& leafs)
{
        leafs.insert(this);
}

// pt_root_t
////////////////////////////////////////////////////////////////////////////////

pt_root_t::pt_root_t(pt_node_t* left,
                     pt_node_t* right,
                     const std::string name,
                     double d)
        : pt_node_t(d, left, right, name, true),
          outgroup_x(-1), outgroup_name("")
{

        id_t leaf_id = 0;
        id_t node_id = n_leafs;
        set_id(this, leaf_id, node_id);

        leaf_map = leaf_map_t(n_leafs, (pt_leaf_t*)NULL);
        node_map = node_map_t(n_nodes, (pt_node_t*)NULL);
        create_leaf_map(this);
        create_node_map(this);
}

pt_root_t::pt_root_t(const pt_root_t& root)
        : pt_node_t(root), leaf_map(root.leaf_map)
{
        // recreate leaf map since we clone all leafs
        create_leaf_map(this);
        create_node_map(this);
}

pt_root_t*
pt_root_t::clone() const
{
        return new pt_root_t(*this);
}

pt_root_t::leafs_t
pt_root_t::get_leafs()
{
        leafs_t leafs;
        get_leafs(leafs);
        return leafs;
}

pt_root_t::nodes_t
pt_root_t::get_nodes()
{
        nodes_t nodes;
        get_nodes(nodes);
        return nodes;
}

pt_root_t::id_t
pt_root_t::get_id(const std::string& taxon) const
{
        for (node_map_t::const_iterator it = node_map.begin(); it != node_map.end(); it++) {
                if ((*it)->name == taxon) {
                        return (*it)->id;
                }
        }
        return -1;
}

bool
pt_root_t::outgroup() const
{
        return d != -HUGE_VAL;
}

void
pt_root_t::create_leaf_map(pt_node_t* node)
{
        if (node->leaf()) {
                pt_leaf_t* leaf = static_cast<pt_leaf_t*>(node);
                leaf_map[leaf->id] = leaf;
        }
        else {
                create_leaf_map(node->left );
                create_leaf_map(node->right);
        }
}

void
pt_root_t::create_node_map(pt_node_t* node)
{
        if (node->leaf()) {
                node_map[node->id] = node;
        }
        else {
                create_node_map(node->left );
                create_node_map(node->right);
        }
}

void
pt_root_t::set_id(pt_node_t* node, pt_node_t::id_t& leaf_id, pt_node_t::id_t& node_id)
{
        if (node->leaf()) {
                node->id = leaf_id++;
        }
        else {
                node->id = node_id++;
                set_id(node->left , leaf_id, node_id);
                set_id(node->right, leaf_id, node_id);
        }
}

// ostream
////////////////////////////////////////////////////////////////////////////////

static const short TAB_WIDTH = 2;

static void
indent(std::ostream& o, size_t nesting)
{
        for(size_t i = 0; i < nesting*TAB_WIDTH; i++) {
                o << " ";
        }
}

static ostream&
print_phylotree(ostream& o, const pt_node_t* node, size_t nesting)
{
        if (node->root()) {
                if (node->id != -1) {
                        o << "(root:"
                          << node->id << " ";
                }
                else {
                        o << "(root ";
                }
                if (!node->leaf()) {
                        print_phylotree(o, node->left,  nesting+1);
                        print_phylotree(o, node->right, nesting+1);
                }
                o << ")"
                  << std::endl;
        }
        else if (node->leaf()) {
                o << std::endl;
                indent(o, nesting);
                if (node->id != -1) {
                        o << "(" << node->name
                          << ":" << node->id
                          << " " << node->d;
                }
                else {
                        o << "(" << node->name
                          << " " << node->d;
                }
                if (static_cast<const pt_leaf_t*>(node)->x != -1) {
                        o << " "
                          << static_cast<const pt_leaf_t*>(node)->x;
                }
                o << ")";
        }
        else {
                o << std::endl;
                indent(o, nesting);
                if (node->id != -1) {
                        o << "(node:"
                          << node->id << " "
                          << node->d;
                }
                else {
                        o << "(node "
                          << node->d;
                }
                print_phylotree(o, node->left,  nesting+1);
                print_phylotree(o, node->right, nesting+1);
                o << ")";
        }
        
        return o;
}

static ostream&
print_newick(ostream& o, const pt_node_t* node, size_t nesting)
{
        if (node->root()) {
                const pt_root_t* root = static_cast<const pt_root_t*>(node);
                o << "(";
                print_newick(o, node->left, nesting+1);
                o << ",";
                print_newick(o, node->right, nesting+1);
                if (root->outgroup()) {
                        o << ","
                          << root->outgroup_name
                          << ":"
                          << root->d;
                }
                o << ");";
        }
        else if (node->leaf()) {
                o << node->name
                  << ":"
                  << node->d;
        }
        else {
                o << "(";
                print_newick(o, node->left, nesting+1);
                o << ",";
                print_newick(o, node->right, nesting+1);
                o << "):"
                  << node->d;
        }
        return o;
}

static std::ostream&
print_newick(std::ostream& o, const pt_node_t* node)
{
        return print_newick(o, node, 0);
}


newick_format::newick_format(const pt_root_t* tree)
        : tree(tree)
{ }

std::ostream&
newick_format::operator()(std::ostream& o) const
{
        return print_newick(o, tree);
}

ostream& operator<< (ostream& o, const pt_node_t* node)
{
        print_phylotree(o, node, 0);
        return o;
}

ostream& operator<< (ostream& o, const newick_format& nf)
{
        return nf(o);
}
