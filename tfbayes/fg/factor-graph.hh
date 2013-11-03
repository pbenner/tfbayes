/* Copyright (C) 2013 Philipp Benner
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

#ifndef FG_FACTOR_GRAPH_HH
#define FG_FACTOR_GRAPH_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <node-types.hh>

class factor_graph_t {
public:
        factor_graph_t(const std::vector<variable_node_i*> variable_nodes,
                       const std::vector<  factor_node_i*> factor_nodes) :
                _variable_nodes(variable_nodes),
                _factor_nodes  (factor_nodes)
                { }

protected:
        std::vector<variable_node_i*> _variable_nodes;
        std::vector<  factor_node_i*> _factor_nodes;
};

#endif /* FG_FACTOR_GRAPH_HH */
