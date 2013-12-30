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
#include <limits>

#include <phylotree.hh>

using boost::optional;

using namespace std;

// pt_node_t
////////////////////////////////////////////////////////////////////////////////

pt_node_t::pt_node_t(double d,
                     pt_node_t* left,
                     pt_node_t* right,
                     const string& name,
                     id_t id)
        : d        (d),
          name     (name),
          id       (id),
          _left    (left),
          _right   (right),
          _ancestor(NULL)
{
        assert ((left && right) || (!left && !right));

        if (leaf()) {
                n_leaves = 1;
                n_nodes  = 1;
        }
        else {
                // update leaf counts and node counts
                n_leaves = left->n_leaves + right->n_leaves;
                n_nodes  = left->n_nodes  + right->n_nodes + 1;
                // record this node as the ancestor
                this->left ()._ancestor = this;
                this->right()._ancestor = this;
        }
}

pt_node_t::pt_node_t(const pt_node_t& node)
        : d        (node.d),
          name     (node.name),
          n_leaves (node.n_leaves),
          n_nodes  (node.n_nodes),
          id       (node.id),
          _left    (NULL),
          _right   (NULL),
          _ancestor(NULL)
{
        if (!node.leaf()) {
                _left  = node.left() .clone();
                _right = node.right().clone();
                left ()._ancestor = this;
                right()._ancestor = this;
        }
}

pt_node_t::~pt_node_t()
{
        if (!leaf()) {
                delete(_left );
                delete(_right);
        }
}

pt_node_t*
pt_node_t::clone() const
{
        return new pt_node_t(*this);
}

void swap(pt_node_t& first, pt_node_t& second)
{
        swap(first.d,         second.d);
        swap(first.name,      second.name);
        swap(first.n_leaves,  second.n_leaves);
        swap(first.n_nodes,   second.n_nodes);
        swap(first.id,        second.id);
        swap(first._left,     second._left);
        swap(first._right,    second._right);
        swap(first._ancestor, second._ancestor);
        // relink ancestors
        if (!first.leaf()) {
                first._left ->_ancestor = &first;
                first._right->_ancestor = &first;
        }
        if (!second.leaf()) {
                second._left ->_ancestor = &second;
                second._right->_ancestor = &second;
        }
}

pt_node_t&
pt_node_t::operator=(const pt_node_t& pt_node)
{
        pt_node_t tmp(pt_node);
        swap(*this, tmp);
        return *this;
}

bool
pt_node_t::leaf() const
{
        return _left == NULL && _right == NULL;
}

bool
pt_node_t::root() const
{
        return _ancestor == NULL;
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
                left ().scale(c);
                right().scale(c);
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

        left ().create_mappings(leaf_map, leaves, node_map, nodes);
        right().create_mappings(leaf_map, leaves, node_map, nodes);
}

void
pt_node_t::set_id(pt_node_t::id_t& leaf_id, pt_node_t::id_t& node_id)
{
        left ().set_id(leaf_id, node_id);
        right().set_id(leaf_id, node_id);
        id = node_id++;
}

void
pt_node_t::set_id(const pt_root_t& tree, pt_node_t::id_t& node_id)
{
        left ().set_id(tree, node_id);
        right().set_id(tree, node_id);
        id = node_id++;
}

void
pt_node_t::move_a()
{
        pt_node_t* tmp;

        if (root() || leaf()) {
                return;
        }
        if (ancestor()._left == this) {
                // relink ancestors
                left()._ancestor  = _ancestor;
                // relink node
                tmp               = _left;
                _left             = ancestor()._right;
                ancestor()._right = tmp;
                // relink ancestors
                left()._ancestor  = this;
        }
        else {
                // relink ancestors
                right()._ancestor = _ancestor;
                tmp               = _right;
                _right            = ancestor()._left;
                ancestor()._left  = tmp;
                // relink ancestors
                right()._ancestor = this;
        }
}

void
pt_node_t::move_b()
{
        pt_node_t* tmp;

        if (root() || leaf()) {
                return;
        }
        if (ancestor()._left == this) {
                // relink ancestors
                _right->_ancestor = _ancestor;
                // relink node
                tmp               = _right;
                _right            = _ancestor->_right;
                _ancestor->_right = tmp;
                // relink ancestors
                _right->_ancestor = this;
        }
        else {
                // relink ancestors
                _left->_ancestor  = _ancestor;
                // relink node
                tmp               = _left;
                _left             = _ancestor->_left;
                _ancestor->_left  = tmp;
                // relink ancestors
                _left->_ancestor  = this;
        }
}

void
pt_node_t::move(size_t which)
{
        switch (which) {
        case 0: move_a(); break;
        case 1: move_b(); break;
        default: break;
        }
}

// pt_leaf_t
////////////////////////////////////////////////////////////////////////////////

pt_leaf_t::pt_leaf_t()
        : pt_node_t()
{ }

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
pt_leaf_t::set_id(const pt_root_t& tree, pt_node_t::id_t& node_id)
{
        assert(name != "");
        id = tree.get_leaf_id(name);
        assert(id != -1);
}

// pt_root_t
////////////////////////////////////////////////////////////////////////////////

pt_root_t::pt_root_t(pt_node_t* left,
                     pt_node_t* right,
                     pt_leaf_t* outgroup,
                     const std::string name,
                     optional<const pt_root_t&> tree)
        : pt_node_t(-numeric_limits<double>::infinity(), left, right, name),
          _outgroup(outgroup)
{
        assert(left  != NULL);
        assert(right != NULL);
        // count outgroup as leaf if present and set the ancestor of
        // the outgroup
        if (outgroup) {
                n_leaves++; n_nodes++;
                this->outgroup()->_ancestor = this;
        }
        // copy leaf ids from second tree
        if (tree) {
                assert (n_leaves == tree->n_leaves);
                id_t node_id = n_leaves;
                set_id(*tree, node_id);
        }
        // set leaf/node ids
        else {
                id_t leaf_id = 0;
                id_t node_id = n_leaves;
                set_id(leaf_id, node_id);
        }

        leaves = leaves_t(n_leaves, (pt_leaf_t*)NULL);
        nodes  = nodes_t (n_nodes,  (pt_node_t*)NULL);
        create_mappings();
}

pt_root_t::pt_root_t(
        pt_node_t& node,
        pt_leaf_t* outgroup,
        optional<const pt_root_t&> tree)
{
        new (this) pt_root_t(node.left().clone(), node.right().clone(),
                             outgroup, node.name, tree);
        delete(&node);
}

pt_root_t::pt_root_t(const pt_root_t& root)
        : pt_node_t(root),
          // copy leaves and nodes since those are vectors
          leaves(root.leaves),
          nodes(root.nodes),
          _outgroup(NULL)
{
        if (root.outgroup()) {
                this->_outgroup = root.outgroup()->clone();
                this-> outgroup()->_ancestor = this;
        }
        // recreate all mappings since we clone all nodes
        create_mappings();
}

pt_root_t::~pt_root_t()
{
        if (_outgroup != NULL) {
                delete(_outgroup);
        }
}

pt_root_t*
pt_root_t::clone() const
{
        return new pt_root_t(*this);
}

void swap(pt_root_t& first, pt_root_t&second)
{
        swap(static_cast<pt_node_t&>(first),
             static_cast<pt_node_t&>(second));
        swap(first.leaf_map,  second.leaf_map);
        swap(first.node_map,  second.node_map);
        swap(first.leaves,    second.leaves);
        swap(first.nodes,     second.nodes);
        swap(first._outgroup, second._outgroup);
        // update maps
        if (first .name != "") first .node_map[first .name] = &first;
        if (second.name != "") second.node_map[second.name] = &second;
        first .nodes[first .id] = &first;
        second.nodes[second.id] = &second;
}

pt_root_t&
pt_root_t::operator=(const pt_root_t& pt_root)
{
        pt_root_t tmp(pt_root);
        swap(*this, tmp);
        return *this;
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

optional<pt_leaf_t&>
pt_root_t::operator()(const std::string& taxon)
{
        leaf_map_t::iterator it = leaf_map.find(taxon);
        if (it != leaf_map.end()) {
                return *it->second;
        }
        return optional<pt_leaf_t&>();
}

optional<const pt_leaf_t&>
pt_root_t::operator()(const std::string& taxon) const
{
        return operator()(taxon);
}

optional<pt_leaf_t&>
pt_root_t::operator()(id_t id)
{
        return *leaves[id];
}

optional<const pt_leaf_t&>
pt_root_t::operator()(id_t id) const
{
        return operator()(id);
}

void
pt_root_t::create_mappings(leaf_map_t& leaf_map, leaves_t& leaves,
                           node_map_t& node_map, nodes_t& nodes)
{
        pt_node_t::create_mappings(leaf_map, leaves, node_map, nodes);

        // also add the outgroup if available
        if (outgroup()) {
                outgroup()->create_mappings(leaf_map, leaves, node_map, nodes);
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
        if (outgroup()) {
                outgroup()->id = leaf_id++;
        }
        pt_node_t::set_id(leaf_id, node_id);
}

void
pt_root_t::set_id(const pt_root_t& tree, pt_node_t::id_t& node_id)
{
        // check for the outgroup and if it is available make sure
        // that it gets the first id
        if (outgroup()) {
                assert(outgroup()->name != "");
                outgroup()->id = tree.get_leaf_id(outgroup()->name);
                assert(outgroup()->id != -1);
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
print_phylotree(ostream& o, const pt_node_t& node, size_t nesting)
{
        if (node.root()) {
                if (node.id != -1) {
                        o << "(root:"
                          << node.id << " ";
                }
                else {
                        o << "(root ";
                }
                if (!node.leaf()) {
                        print_phylotree(o, node.left (), nesting+1);
                        print_phylotree(o, node.right(), nesting+1);
                }
                o << ")"
                  << std::endl;
        }
        else if (node.leaf()) {
                o << std::endl;
                indent(o, nesting);
                if (node.id != -1) {
                        o << "(" << node.name
                          << ":" << node.id
                          << " " << node.d;
                }
                else {
                        o << "(" << node.name
                          << " " << node.d;
                }
                o << ")";
        }
        else {
                o << std::endl;
                indent(o, nesting);
                if (node.id != -1) {
                        o << "(node:"
                          << node.id << " "
                          << node.d;
                }
                else {
                        o << "(node "
                          << node.d;
                }
                print_phylotree(o, node.left (), nesting+1);
                print_phylotree(o, node.right(), nesting+1);
                o << ")";
        }
        
        return o;
}

static ostream&
print_newick(ostream& o, const pt_node_t& node, size_t nesting)
{
        if (node.root()) {
                const pt_root_t& root = static_cast<const pt_root_t&>(node);
                o << "(";
                print_newick(o, node.left (), nesting+1);
                o << ",";
                print_newick(o, node.right(), nesting+1);
                if (root.outgroup()) {
                        o << ","
                          << root.outgroup()->name
                          << ":"
                          << fixed << root.outgroup()->d;
                }
                o << ");";
        }
        else if (node.leaf()) {
                o << node.name
                  << ":"
                  << fixed << node.d;
        }
        else {
                if (node.d != 0.0) {
                        o << "(";
                        print_newick(o, node.left (), nesting+1);
                        o << ",";
                        print_newick(o, node.right(), nesting+1);
                        o << "):"
                          << fixed << node.d;
                }
                else {
                        print_newick(o, node.left (), nesting+1);
                        o << ",";
                        print_newick(o, node.right(), nesting+1);
                }
        }
        return o;
}

static std::ostream&
print_newick(std::ostream& o, const pt_node_t& node)
{
        return print_newick(o, node, 0);
}


newick_format::newick_format(const pt_root_t& tree)
        : tree(&tree)
{ }

std::ostream&
newick_format::operator()(std::ostream& o) const
{
        return print_newick(o, *tree);
}

ostream& operator<< (ostream& o, const pt_node_t& node)
{
        return print_phylotree(o, node, 0);
}

ostream& operator<< (ostream& o, const newick_format& nf)
{
        return nf(o);
}
