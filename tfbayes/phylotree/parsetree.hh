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

#ifndef __TFBAYES_PHYLOTREE_PARSETREE_HH__
#define __TFBAYES_PHYLOTREE_PARSETREE_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <cstdarg>
#include <ostream>

#include <tfbayes/phylotree/phylotree.hh>

typedef enum {
        /* tree list */
        TREE_LIST_N, TREE_N,

        /* nodes and leaf */
        NODE_LIST_N, NODE_N, LEAF_N,

        /* Species name */
        NAME_N,

        /* Distance to ancestor */
        DISTANCE_N,
} nodetype_t;

class pt_parsetree_t {
public:
         pt_parsetree_t(nodetype_t type, size_t n_children, void *data, ...);
        ~pt_parsetree_t();

        std::list<pt_root_t> convert(
                size_t drop, size_t skip,
                boost::optional<const pt_root_t&> ref_tree = boost::optional<const pt_root_t&>()) const;
        pt_root_t  convert(
                const std::list<pt_root_t>& tree_list,
                boost::optional<const pt_root_t&> ref_tree) const;
        pt_node_t* convert() const;

        void * const data;
        const nodetype_t type;
        const size_t n_children;
        pt_parsetree_t **children;
};

std::list<pt_root_t> parse_tree_list(FILE * file = NULL, size_t drop = 0, size_t skip = 1,
                                     boost::optional<const pt_root_t&> ref_tree = boost::optional<const pt_root_t&>());
std::list<pt_root_t> parse_tree_list(const std::string& filename, size_t drop = 0, size_t skip = 1,
                                     boost::optional<const pt_root_t&> ref_tree = boost::optional<const pt_root_t&>());

std::ostream& operator<< (std::ostream& o, pt_parsetree_t* const tree);

#endif /* __TFBAYES_PHYLOTREE_PARSETREE_HH__ */
