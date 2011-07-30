/* Copyright (C) 2011 Philipp Benner
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

#ifndef DPM_TFBS_GRAPH_HH
#define DPM_TFBS_GRAPH_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <map>

#include <index.hh>

typedef struct _edge_t {
        bool operator==(const _edge_t& edge) const {
                if ((index1 == edge.index1 && index2 == edge.index2) ||
                    (index1 == edge.index2 && index2 == edge.index1)) {
                        return true;
                }
                return false;

        }
        const index_t& index1;
        const index_t& index2;
} edge_t;

class Graph {
public:
        Graph() {}

protected:
        std::map<edge_t, size_t> _edges;
};

#endif /* DPM_TFBS_GRAPH_HH */
