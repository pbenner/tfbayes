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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <phylotree.hh>

pt_node_t::pt_node_t(short x, double d, pt_node_t* left, pt_node_t* right) {
        this->x = x;
        this->d = d;
        this->left  = left;
        this->right = right;
}

bool
pt_node_t::leaf() const {
        return left == NULL && right == NULL;
}

bool
pt_node_t::root() const {
        return d == 0.0;
}
