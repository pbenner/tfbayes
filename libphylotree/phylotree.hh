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

#ifndef PHYLOTREE_HH
#define PHYLOTREE_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <math.h>
#include <cstddef>
#include <set>
#include <string>
#include <ostream>

#include <clonable.hh>

#include <boost/unordered_map.hpp>

typedef short nucleotide_t;

class pt_node_t : public clonable {
public:
        pt_node_t(nucleotide_t x = -1, double d = 0.0,
                  pt_node_t* left  = NULL,
                  pt_node_t* right = NULL,
                  const std::string& name = "",
                  const ssize_t id = -1)
                : x(x), d(d), left(left), right(right), name(name), id(id)
                { }
        pt_node_t(const pt_node_t& node)
                : x(node.x), d(node.d), left(NULL), right(NULL), name(node.name), id(node.id) {

                if (!node.leaf()) {
                        left  = node.left ->clone();
                        right = node.right->clone();
                }
        }
        pt_node_t* clone() const {
                return new pt_node_t(*this);
        }

        /* typedef for a set of nodes */
        typedef std::set<pt_node_t*> nodes_t;
        typedef ssize_t id_t;

        void destroy() {
                if (!leaf()) {
                        left ->destroy();
                        right->destroy();
                }
                delete(this);
        }

        bool leaf() const { return left == NULL && right == NULL; }
        bool root() const { return d == -HUGE_VAL; }

        void init(const nucleotide_t& nucleotide) {
                this->x = nucleotide;
                if (!leaf()) {
                        left ->init(nucleotide);
                        right->init(nucleotide);
                }
        }

        nodes_t get_nodes() {
                nodes_t nodes;
                get_nodes(nodes);
                return nodes;
        }

        double mutation_probability() const {
                return 1.0-exp(-d);
        }
        /* scale distances in the tree by a constant factor */
        void scale(double c) {
                d *= c;
                if (!leaf()) {
                        left->scale(c);
                        right->scale(c);
                }
        }

        std::ostream& print_distances(std::ostream& o) {
                if (root()) { 
                        o << "{";
                        left ->print_distances(o);
                        o  << ", ";
                        right->print_distances(o);
                        o << "}";
                }
                else {
                        o << d;
                        if (!leaf()) {
                                o  << ", ";
                                left ->print_distances(o);
                                o  << ", ";
                                right->print_distances(o);
                        }
                }
                return o;
        }
        std::ostream& print_phylotree(std::ostream& o, size_t nesting) const {
                if (root()) {
                        o << "(root";
                        if (id != -1) {
                                o << ":" << id;
                        }
                        if (!leaf()) {
                                left ->print_phylotree(o, nesting+1);
                                right->print_phylotree(o, nesting+1);
                        }
                        o << ")"
                          << std::endl;
                }
                else if (leaf()) {
                        o << std::endl;
                        indent(o, nesting);
                        if (id != -1) {
                                o << "(" << name
                                  << ":" << id
                                  << " " << d;
                        }
                        else {
                                o << "(" << name
                                  << " " << d;
                        }
                        if (x != -1) {
                                o << " "
                                  << x;
                        }
                        o << ")";
                }
                else {
                        o << std::endl;
                        indent(o, nesting);
                        if (id != -1) {
                                o << "(node:"
                                  << id << " "
                                  << d;
                        }
                        else {
                                o << "(node "
                                  << d;
                        }
                        left ->print_phylotree(o, nesting+1);
                        right->print_phylotree(o, nesting+1);
                        o << ")";
                }
                return o;
        }
        std::ostream& print_phylotree(std::ostream& o) const {
                return print_phylotree(o, 0);
        }
        std::ostream& print_newick(std::ostream& o, size_t nesting) const {
                if (root()) {
                        o << "(";
                        left ->print_newick(o, nesting+1);
                        o << ",";
                        right->print_newick(o, nesting+1);
                        o << ");"
                          << std::endl;
                }
                else if (leaf()) {
                        o << std::endl;
                        indent(o, nesting);
                        o << name
                          << ":"
                          << d;
                }
                else {
                        o << "(";
                        left ->print_newick(o, nesting+1);
                        o << ",";
                        right->print_newick(o, nesting+1);
                        o << "):"
                          << d;
                }
                return o;
        }
        std::ostream& print_newick(std::ostream& o) const {
                return print_newick(o, 0);
        }
        std::ostream& print(std::ostream& o, bool newick = false) const {
                if (newick) {
                        return print_newick(o);
                }
                else {
                        return print_phylotree(o);
                }
                return o;
        }
        id_t get_id(const std::string& taxon) {
                id_t tmp;
                if (leaf()) {
                        if (name == taxon) {
                                return id;
                        }
                        else {
                                return -1;
                        }
                }
                else {
                        if ((tmp = left->get_id(taxon)) != -1) {
                                return tmp;
                        }
                        else {
                                return right->get_id(taxon);
                        }
                }
        }

        /* coded nucleotide */
        nucleotide_t x;
        /* distance to ancestor */
        double d;
        /* left child */
        pt_node_t* left;
        /* right child */
        pt_node_t* right;
        /* name of the node */
        const std::string name;
        /* identifier */
        id_t id;

private:
        void get_nodes(nodes_t& nodes) {
                if (!root()) {
                        nodes.insert(this);
                }
                if (!leaf()) {
                        left ->get_nodes(nodes);
                        right->get_nodes(nodes);
                }
        }
        static const short TAB_WIDTH = 2;
        static
        void indent(std::ostream& o, size_t nesting) {
                for(size_t i = 0; i < nesting*TAB_WIDTH; i++) {
                        o << " ";
                }
        }
};

#include <vector>

class pt_root_t : public pt_node_t {
public:
        /* typedefs */
        typedef std::vector<pt_node_t*> node_map_t;

        pt_root_t(short x = -1,
                  pt_node_t* left  = NULL,
                  pt_node_t* right = NULL,
                  const std::string name = "")
                : pt_node_t(x, -HUGE_VAL, left, right, name) {

                pt_node_t::id_t n = set_id(this, 0)+1;
                node_map = node_map_t(n, (pt_node_t*)NULL);
                create_map(this);
        }
        pt_root_t(const pt_root_t& root)
                : pt_node_t(root), node_map(root.node_map) {
                create_map(this);
        }
        pt_root_t* clone() const {
                return new pt_root_t(*this);
        }
        pt_node_t::id_t size() const {
                return node_map.size();
        }

        node_map_t node_map;

protected:
        void create_map(pt_node_t* node) {
                node_map[node->id] = node;
                if (!node->leaf()) {
                        create_map(node->left );
                        create_map(node->right);
                }
        }
        pt_node_t::id_t set_id(pt_node_t* node, pt_node_t::id_t id) {
                node->id = id;
                if (!node->leaf()) {
                        id = set_id(node->left , id+1);
                        id = set_id(node->right, id+1);
                }
                return id;
        }
};

class pt_leaf_t : public pt_node_t {
public:
        pt_leaf_t(short x, double d, const std::string name = "")
                : pt_node_t(x, d, NULL, NULL, name) {
        }
};

std::ostream& operator<< (std::ostream& o, pt_node_t* const node);

#endif /* PHYLOTREE_HH */
