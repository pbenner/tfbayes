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
                     id_t id)
        : d(d), left(left), right(right), ancestor(NULL),
          name(name), id(id)
{
        assert ((left == NULL && right == NULL) ||
                (left != NULL && right != NULL));

        if (leaf()) {
                n_leaves = 1;
                n_nodes  = 1;
        }
        else {
                // update leaf counts and node counts
                n_leaves = left->n_leaves + right->n_leaves;
                n_nodes  = left->n_nodes  + right->n_nodes + 1;
                // record this node as the ancestor
                left ->ancestor = this;
                right->ancestor = this;
        }
}

pt_node_t::pt_node_t(const pt_node_t& node)
        : d(node.d),
          left(NULL), right(NULL),
          ancestor(NULL),
          name(node.name),
          n_leaves(node.n_leaves),
          n_nodes(node.n_nodes),
          id(node.id)
{
        
        if (!node.leaf()) {
                left  = node.left ->clone();
                right = node.right->clone();
                left ->ancestor = this;
                right->ancestor = this;
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

pt_root_t*
pt_node_t::convert_to_root(pt_leaf_t* outgroup,
                           const pt_root_t* tree)
{
        pt_root_t* pt_root = new pt_root_t(left, right, outgroup, name, tree);
        delete(this);
        return pt_root;
}

bool
pt_node_t::leaf() const
{
        return left == NULL && right == NULL;
}

bool
pt_node_t::root() const
{
        return ancestor == NULL;
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
pt_node_t::create_mappings(leaf_map_t& leaf_map, leaves_t& leaves,
                           node_map_t& node_map, nodes_t& nodes)
{
        if (name != "") {
                node_map[name] = this;
        }
        nodes[id] = this;

        left ->create_mappings(leaf_map, leaves, node_map, nodes);
        right->create_mappings(leaf_map, leaves, node_map, nodes);
}

void
pt_node_t::set_id(pt_node_t::id_t& leaf_id, pt_node_t::id_t& node_id)
{
        id = node_id++;
        left ->set_id(leaf_id, node_id);
        right->set_id(leaf_id, node_id);
}

void
pt_node_t::set_id(const pt_root_t* tree, pt_node_t::id_t& node_id)
{
        id = node_id++;
        left ->set_id(tree, node_id);
        right->set_id(tree, node_id);
}

void
pt_node_t::move_a()
{
        pt_node_t* tmp;

        if (root() || leaf()) {
                return;
        }
        if (ancestor->left == this) {
                // relink ancestors
                left->ancestor  = ancestor;
                // relink node
                tmp             = left;
                left            = ancestor->right;
                ancestor->right = tmp;
                // relink ancestors
                left->ancestor  = this;
        }
        else {
                // relink ancestors
                right->ancestor = ancestor;
                tmp             = right;
                right           = ancestor->left;
                ancestor->left  = tmp;
                // relink ancestors
                right->ancestor = this;
        }
}

void
pt_node_t::move_b()
{
        pt_node_t* tmp;

        if (root() || leaf()) {
                return;
        }
        if (ancestor->left == this) {
                // relink ancestors
                right->ancestor = ancestor;
                // relink node
                tmp             = right;
                right           = ancestor->right;
                ancestor->right = tmp;
                // relink ancestors
                right->ancestor = this;
        }
        else {
                // relink ancestors
                left->ancestor  = ancestor;
                // relink node
                tmp             = left;
                left            = ancestor->left;
                ancestor->left  = tmp;
                // relink ancestors
                left->ancestor  = this;
        }
}

// pt_leaf_t
////////////////////////////////////////////////////////////////////////////////

pt_leaf_t::pt_leaf_t(double d, const std::string name)
        : pt_node_t(d, NULL, NULL, name)
{ }

pt_leaf_t::pt_leaf_t(const pt_leaf_t& leaf)
        : pt_node_t(leaf)
{ }

pt_leaf_t*
pt_leaf_t::clone() const
{
        return new pt_leaf_t(*this);
}

void
pt_leaf_t::create_mappings(leaf_map_t& leaf_map, leaves_t& leaves,
                           node_map_t& node_map, nodes_t& nodes)
{
        if (name != "") {
                leaf_map[name] = this;
                node_map[name] = this;
        }
        leaves[id] = this;
        nodes [id] = this;
}

void
pt_leaf_t::set_id(pt_node_t::id_t& leaf_id, pt_node_t::id_t& node_id)
{
        id = leaf_id++;
}

void
pt_leaf_t::set_id(const pt_root_t* tree, pt_node_t::id_t& node_id)
{
        assert(name != "");
        id = tree->get_leaf_id(name);
        assert(id != -1);
}

// pt_root_t
////////////////////////////////////////////////////////////////////////////////

pt_root_t::pt_root_t(pt_node_t* left,
                     pt_node_t* right,
                     pt_leaf_t* outgroup,
                     const std::string name,
                     const pt_root_t* tree)
        : pt_node_t(-HUGE_VAL, left, right, name),
          outgroup(outgroup)
{
        // count outgroup as leaf if present and set the ancestor of
        // the outgroup
        if (has_outgroup()) {
                n_leaves++; n_nodes++;
                outgroup->ancestor = this;
        }
        // set leaf/node ids
        if (tree == NULL) {
                id_t leaf_id = 0;
                id_t node_id = n_leaves;
                set_id(leaf_id, node_id);
        }
        // copy leaf ids from second tree
        else {
                assert (n_leaves == tree->n_leaves);
                id_t node_id = n_leaves;
                set_id(tree, node_id);
        }

        leaves = leaves_t(n_leaves, (pt_leaf_t*)NULL);
        nodes  = nodes_t (n_nodes,  (pt_node_t*)NULL);
        create_mappings();
}

pt_root_t::pt_root_t(const pt_root_t& root)
        : pt_node_t(root),
          // copy leaves and nodes since those are vectors
          leaves(root.leaves),
          nodes(root.nodes),
          outgroup(NULL)
{
        if (root.has_outgroup()) {
                outgroup = root.outgroup->clone();
                outgroup->ancestor = this;
        }
        // recreate all mappings since we clone all nodes
        create_mappings();
}

void
pt_root_t::destroy()
{
        if (has_outgroup()) {
                outgroup->destroy();
        }
        pt_node_t::destroy();
}

pt_root_t*
pt_root_t::clone() const
{
        return new pt_root_t(*this);
}

pt_root_t::id_t
pt_root_t::get_node_id(const std::string& taxon) const
{
        node_map_t::const_iterator it = node_map.find(taxon);
        if (it != node_map.end()) {
                return it->second->id;
        }
        return -1;
}

pt_root_t::id_t
pt_root_t::get_leaf_id(const std::string& taxon) const
{
        leaf_map_t::const_iterator it = leaf_map.find(taxon);
        if (it != leaf_map.end()) {
                return it->second->id;
        }
        return -1;
}

pt_leaf_t*
pt_root_t::operator()(const std::string& taxon)
{
        leaf_map_t::iterator it = leaf_map.find(taxon);
        if (it != leaf_map.end()) {
                return it->second;
        }
        return NULL;
}

const pt_leaf_t*
pt_root_t::operator()(const std::string& taxon) const
{
        return operator()(taxon);
}

pt_leaf_t*
pt_root_t::operator()(id_t id)
{
        return leaves[id];
}

const pt_leaf_t*
pt_root_t::operator()(id_t id) const
{
        return operator()(id);
}

bool
pt_root_t::has_outgroup() const
{
        return outgroup != NULL;
}

void
pt_root_t::create_mappings(leaf_map_t& leaf_map, leaves_t& leaves,
                           node_map_t& node_map, nodes_t& nodes)
{
        pt_node_t::create_mappings(leaf_map, leaves, node_map, nodes);

        // also add the outgroup if available
        if (has_outgroup()) {
                outgroup->create_mappings(leaf_map, leaves, node_map, nodes);
        }
}

void
pt_root_t::create_mappings()
{
        create_mappings(leaf_map, leaves, node_map, nodes);
}

void
pt_root_t::set_id(pt_node_t::id_t& leaf_id, pt_node_t::id_t& node_id)
{
        // check for the outgroup and if it is available make sure
        // that it gets the first id
        if (has_outgroup()) {
                outgroup->id = leaf_id++;
        }
        pt_node_t::set_id(leaf_id, node_id);
}

void
pt_root_t::set_id(const pt_root_t* tree, id_t& node_id)
{
        // check for the outgroup and if it is available make sure
        // that it gets the first id
        if (has_outgroup()) {
                assert(outgroup->name != "");
                outgroup->id = tree->get_leaf_id(outgroup->name);
                assert(outgroup->id != -1);
        }
        pt_node_t::set_id(tree, node_id);
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
                if (root->has_outgroup()) {
                        o << ","
                          << root->outgroup->name
                          << ":"
                          << fixed << root->outgroup->d;
                }
                o << ");";
        }
        else if (node->leaf()) {
                o << node->name
                  << ":"
                  << fixed << node->d;
        }
        else {
                if (node->d != 0.0) {
                        o << "(";
                        print_newick(o, node->left, nesting+1);
                        o << ",";
                        print_newick(o, node->right, nesting+1);
                        o << "):"
                          << fixed << node->d;
                }
                else {
                        print_newick(o, node->left, nesting+1);
                        o << ",";
                        print_newick(o, node->right, nesting+1);
                }
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
