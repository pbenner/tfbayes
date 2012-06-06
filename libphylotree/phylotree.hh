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
                  const std::string& name = "")
                : x(x), d(d), left(left), right(right), name(name)
                { }

        pt_node_t* clone() const {
                pt_node_t* pt_node = new pt_node_t(*this);
                if (!leaf()) {
                        pt_node->left  = left ->clone();
                        pt_node->right = right->clone();
                }
                return pt_node;
        }

        /* typedef for a set of nodes */
        typedef std::set<pt_node_t*> nodes_t;

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
                        left ->print_phylotree(o, nesting+1);
                        right->print_phylotree(o, nesting+1);
                        o << ")"
                          << std::endl;
                }
                else if (leaf()) {
                        o << std::endl;
                        indent(o, nesting);
                        o << "("
                          << name
                          << " "
                          << d;
                        if (x != -1) {
                                o << " "
                                  << x;
                        }
                        o << ")";
                }
                else {
                        o << std::endl;
                        indent(o, nesting);
                        o << "(node "
                          << d;
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

class pt_root_t : public pt_node_t {
public:
        pt_root_t(short x = -1,
                  pt_node_t* left  = NULL,
                  pt_node_t* right = NULL,
                  const std::string name = "")
                : pt_node_t(x, -HUGE_VAL, left, right, name) {

                create_map();
        }
        pt_root_t* clone() const {
                pt_root_t* pt_root = new pt_root_t(*this);
                if (!leaf()) {
                        pt_root->left  = left ->clone();
                        pt_root->right = right->clone();
                }
                pt_root->create_map();

                return pt_root;
        }

        /* typedefs */
        typedef boost::unordered_map<const std::string, pt_node_t*> map_t;

        /* methods */
        void create_map() {
                create_map(this);
        }
        void create_map(pt_node_t* node) {
                if (node->leaf() && node->name != "") {
                        map[node->name] = node;
                }
                else {
                        create_map(node->left );
                        create_map(node->right);
                }
        }

        map_t map;
};

class pt_leaf_t : public pt_node_t {
public:
        pt_leaf_t(short x, double d, const std::string name = "")
                : pt_node_t(x, d, NULL, NULL, name) {
        }
};

std::ostream& operator<< (std::ostream& o, pt_node_t* const node);

#endif /* PHYLOTREE_HH */
