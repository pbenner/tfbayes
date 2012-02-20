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
#include <string>

typedef short nucleotide_t;

class pt_node_t {
public:
        pt_node_t(nucleotide_t x = -1, double d = 0.0,
                  pt_node_t* left  = NULL,
                  pt_node_t* right = NULL,
                  const std::string& name = "")
                : x(x), d(d), left(left), right(right), name(name)
                { }

        bool leaf() const { return left == NULL && right == NULL; }
        bool root() const { return d == 0.0; }

        double mutation_probability() const {
                return 1.0-exp(-d);
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
};

class pt_root_t : public pt_node_t {
public:
        pt_root_t(short x = -1,
                  pt_node_t* left  = NULL,
                  pt_node_t* right = NULL,
                  const std::string name = "")
                : pt_node_t(x, 0.0, left, right, name) { }
};

class pt_leaf_t : public pt_node_t {
public:
        pt_leaf_t(short x, double d, const std::string name = "")
                : pt_node_t(x, d, NULL, NULL, name) { }
};

std::ostream& operator<< (std::ostream& o, pt_node_t* const node);

#endif /* PHYLOTREE_HH */
