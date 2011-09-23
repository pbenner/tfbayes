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

#ifndef GRAPH_HH
#define GRAPH_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <map>
#include <boost/unordered_map.hpp> 

#include <index.hh>
//#include <indexer.hh>
//#include <datatypes.hh>

static inline
const index_t& min(const index_t& index1, const index_t& index2) {
        if ((index1[0] <  index2[0]) ||
            (index1[0] == index2[0] && index1[1] <= index2[1])) {
                return index1;
        }
        else {
                return index2;
        }
}

static inline
const index_t& max(const index_t& index1, const index_t& index2) {
        if ((index1[0] >  index2[0]) ||
            (index1[0] == index2[0] && index1[1] > index2[1])) {
                return index1;
        }
        else {
                return index2;
        }
}

typedef struct _edge_t {
        _edge_t(const index_t& i1, const index_t& i2)
                : index1(min(i1,i2)), index2(max(i1,i2)) { }

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

static inline
size_t hash_value(const edge_t& edge)
{
        boost::hash<size_t> hasher;

        return hasher(edge.index1[0]*edge.index1[1]) + hasher(edge.index2[1]);
}

class Graph {
public:
        Graph() {}

        // typedefs
        ////////////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_map<edge_t, size_t> map_t;

        typedef map_t::iterator iterator;
        typedef map_t::const_iterator const_iterator;

        // iterators
        ////////////////////////////////////////////////////////////////////////////////
        iterator begin() { return _edges.begin(); }
        iterator end()   { return _edges.end();   }

        const_iterator begin() const { return _edges.begin(); }
        const_iterator end()   const { return _edges.end();   }

        // operators
        ////////////////////////////////////////////////////////////////////////////////
        size_t& operator[](const edge_t& edge) {
                return _edges[edge];
        }

        // methods
        ////////////////////////////////////////////////////////////////////////////////
        void insert(const edge_t& edge) {
                _edges[edge]++;
        }
        void insert(const Graph& graph) {
                for (const_iterator it = graph.begin(); it != graph.end(); it++) {
                        _edges[(*it).first] += (*it).second;
                }
        }
        void cleanup(size_t threshold) {
                for (map_t::const_iterator it = _edges.begin();
                     it != _edges.end();) {
                        if ((*it).second <= threshold) {
                                it = _edges.erase(it);
                        }
                        else {
                                it++;
                        }
                }
        }
        void cleanup() {
                for (map_t::const_iterator it = _edges.begin();
                     it != _edges.end();) {
                        it = _edges.erase(it);
                }
        }

protected:
        map_t _edges;
};

#endif /* GRAPH_HH */
