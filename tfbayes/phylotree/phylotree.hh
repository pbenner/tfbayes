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

#ifndef PHYLOTREE_HH
#define PHYLOTREE_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <math.h>
#include <cstddef>
#include <cassert>
#include <set>
#include <string>
#include <ostream>
#include <vector>

#include <boost/unordered_map.hpp>

#include <tfbayes/utility/clonable.hh>

typedef short nucleotide_t;

class pt_leaf_t;
class pt_root_t;

class pt_node_t : public clonable {
public:
        /* typedefs */
        typedef std::vector<pt_leaf_t*> leaf_map_t;
        typedef std::vector<pt_node_t*> node_map_t;
        typedef std::set<pt_leaf_t*> leafs_t;
        typedef std::set<pt_node_t*> nodes_t;
        typedef ssize_t id_t;

        pt_node_t(double d = 0.0,
                  pt_node_t* left  = NULL,
                  pt_node_t* right = NULL,
                  const std::string& name = "",
                  bool is_root = false,
                  id_t id = -1);
        pt_node_t(const pt_node_t& node);

        pt_node_t* clone() const;
        void destroy();

        bool leaf() const;
        bool root() const;

        double mutation_probability() const;
        /* scale distances in the tree by a constant factor */
        void scale(double c);

        /* distance to ancestor */
        double d;
        /* left child */
        pt_node_t* left;
        /* right child */
        pt_node_t* right;
        /* name of the node */
        const std::string name;
        /* is this the root node? */
        bool is_root;
        /* number of leafs in this tree */
        ssize_t n_leafs;
        ssize_t n_nodes;
        /* identifier */
        id_t id;

protected:
        virtual void get_leafs(leafs_t& leafs);
        virtual void get_nodes(nodes_t& nodes);

};

class pt_leaf_t : public pt_node_t {
public:
        pt_leaf_t(short x, double d, const std::string name = "");

        virtual void get_leafs(leafs_t& leafs);

        /* coded nucleotide */
        nucleotide_t x;
};

class pt_root_t : public pt_node_t {
public:
        pt_root_t(pt_node_t* left  = NULL,
                  pt_node_t* right = NULL,
                  const std::string name = "",
                  double d = -HUGE_VAL);
        pt_root_t(const pt_root_t& root);

        pt_root_t* clone() const;

        using pt_node_t::get_nodes;
        using pt_node_t::get_leafs;
        virtual nodes_t get_nodes();
        virtual leafs_t get_leafs();

        id_t get_id(const std::string& taxon) const;
        
        leaf_map_t leaf_map;
        node_map_t node_map;

protected:
        void create_leaf_map(pt_node_t* node);
        void create_node_map(pt_node_t* node);
        void set_id(pt_node_t* node, id_t& leaf_id, id_t& node_id);
};

class newick_format {
public:
        newick_format(const pt_root_t* tree);
        std::ostream& operator()(std::ostream& o) const;
protected:
        const pt_root_t* tree;
};

std::ostream& operator<< (std::ostream& o, const pt_node_t* node);
std::ostream& operator<< (std::ostream& o, const newick_format& nf);

#endif /* PHYLOTREE_HH */
