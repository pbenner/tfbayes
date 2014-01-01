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

#ifndef __TFBAYES_PHYLOTREE_HH__
#define __TFBAYES_PHYLOTREE_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <math.h>
#include <cstddef>
#include <cassert>
#include <string>
#include <ostream>
#include <vector>

#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>

#include <tfbayes/utility/clonable.hh>

class pt_leaf_t;
class pt_root_t;

class pt_node_t : public virtual clonable {
public:
        /* typedefs */
        typedef boost::unordered_map<std::string, pt_leaf_t*> leaf_map_t;
        typedef boost::unordered_map<std::string, pt_node_t*> node_map_t;
        typedef std::vector<pt_leaf_t*> leaves_t;
        typedef std::vector<pt_node_t*> nodes_t;
        typedef ssize_t id_t;

        pt_node_t(double d = 0.0,
                  pt_node_t* left  = NULL,
                  pt_node_t* right = NULL,
                  const std::string& name = "",
                  id_t id = -1);
        explicit pt_node_t(const pt_node_t& node);
        virtual ~pt_node_t();

        pt_node_t* clone() const;

        friend void swap(pt_node_t& first, pt_node_t&second);

        virtual pt_node_t& operator=(const pt_node_t& pt_node);

        bool leaf() const;
        bool root() const;

        const pt_node_t& left () const {
                return *_left;
        }
              pt_node_t& left () {
                return *_left;
        }
        const pt_node_t& right() const {
                return *_right;
        }
              pt_node_t& right() {
                return *_right;
        }
        const pt_node_t& ancestor() const {
                return *_ancestor;
        }
              pt_node_t& ancestor() {
                return *_ancestor;
        }

        double mutation_probability() const;
        /* scale distances in the tree by a constant factor */
        void scale(double c);

        /* moves */
        void move_a();
        void move_b();
        void move(size_t which);

        /* distance to ancestor */
        double d;
        /* name of the node */
        std::string name;
        /* number of leaves in this tree */
        ssize_t n_leaves;
        ssize_t n_nodes;
        /* identifier */
        id_t id;

        friend class pt_root_t;
        friend class pt_leaf_t;

protected:
        virtual void create_mappings(leaf_map_t& leaf_map, leaves_t& leaves,
                                     node_map_t& node_map, nodes_t& nodes);

        virtual void set_id(id_t& leaf_id, id_t& node_id);
        virtual void set_id(const pt_root_t& tree, id_t& node_id);

        /* left child */
        pt_node_t* _left;
        /* right child */
        pt_node_t* _right;
        /* link to ancestor */
        pt_node_t* _ancestor;
};

class pt_leaf_t : public pt_node_t {
public:
        pt_leaf_t();
        pt_leaf_t(double d, const std::string name = "");
        pt_leaf_t(const pt_leaf_t& leaf);

        pt_leaf_t* clone() const;

        using pt_node_t::operator=;

        virtual void create_mappings(leaf_map_t& leaf_map, leaves_t& leaves,
                                     node_map_t& node_map, nodes_t& nodes);

protected:
        virtual void set_id(id_t& leaf_id, id_t& node_id);
        virtual void set_id(const pt_root_t& tree, id_t& node_id);
};

class pt_root_t : public pt_node_t {
public:
        pt_root_t(pt_node_t* left,
                  pt_node_t* right,
                  pt_leaf_t* outgroup = NULL,
                  const std::string name = "",
                  // if root is not null then copy ids from this tree
                  boost::optional<const pt_root_t&> tree = boost::optional<const pt_root_t&>());
        // convert a normal node to a root (and delete the old node)
        pt_root_t(pt_node_t& node,
                  pt_leaf_t* outgroup = NULL,
                  // if root is not null then copy ids from this tree
                  boost::optional<const pt_root_t&> tree = boost::optional<const pt_root_t&>());
        pt_root_t(const pt_root_t& root);
        virtual ~pt_root_t();

        pt_root_t* clone() const;

        friend void swap(pt_root_t& first, pt_root_t&second);

        pt_root_t& operator=(const pt_root_t& pt_root);

        id_t get_node_id(const std::string& taxon) const;
        id_t get_leaf_id(const std::string& taxon) const;

        boost::optional<const pt_leaf_t&> outgroup() const {
                return _outgroup ? *_outgroup : boost::optional<const pt_leaf_t&>();
        }
        boost::optional<      pt_leaf_t&> outgroup() {
                return _outgroup ? *_outgroup : boost::optional<      pt_leaf_t&>();
        }

        boost::optional<      pt_leaf_t&> operator()(const std::string& taxon);
        boost::optional<const pt_leaf_t&> operator()(const std::string& taxon) const;
        boost::optional<      pt_leaf_t&> operator()(id_t id);
        boost::optional<const pt_leaf_t&> operator()(id_t id) const;

        // leaf or node name -> leaf or node
        leaf_map_t leaf_map;
        node_map_t node_map;
        // leaf or node id -> leaf or node
        leaves_t leaves;
        nodes_t nodes;

protected:
        virtual void create_mappings(leaf_map_t& leaf_map, leaves_t& leaves,
                                     node_map_t& node_map, nodes_t& nodes);
        virtual void create_mappings();

        virtual void set_id(id_t& leaf_id, id_t& node_id);
        virtual void set_id(const pt_root_t& tree, id_t& node_id);

        pt_leaf_t* _outgroup;
};

class newick_format {
public:
        newick_format(const pt_root_t& tree);
        std::ostream& operator()(std::ostream& o) const;
protected:
        const pt_root_t* tree;
};

std::ostream& operator<< (std::ostream& o, const pt_node_t& node);
std::ostream& operator<< (std::ostream& o, const newick_format& nf);

#endif /* __TFBAYES_PHYLOTREE_HH__ */
