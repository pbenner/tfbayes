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
#include <tr1/unordered_map>

#include <index.hh>
#include <indexer.hh>
#include <datatypes.hh>

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

namespace std {
namespace tr1 {
        template<>
        class hash<edge_t>
                : public unary_function<edge_t, size_t>
        {
        public:
                size_t operator()(const edge_t& edge) const {
                        return static_cast<size_t>(edge.index1[0] << sizeof(size_t)*4 & edge.index1[1]);
                }
        };
}
}

class TfbsGraph {
public:
        TfbsGraph() {}

        // typedefs
        ////////////////////////////////////////////////////////////////////////////////
        typedef std::tr1::unordered_map<edge_t, size_t> map_t;

        typedef map_t::iterator iterator;
        typedef map_t::const_iterator const_iterator;

        // iterators
        ////////////////////////////////////////////////////////////////////////////////
        iterator begin() { return _edges.begin(); }
        iterator end()   { return _edges.end(); }

        const_iterator begin() const { return _edges.begin(); }
        const_iterator end()   const { return _edges.end(); }

        // methods
        ////////////////////////////////////////////////////////////////////////////////
        void insert(const edge_t& edge) {
                _edges[edge]++;
        }
        void update(const Indexer& indexer,
                    const sequence_data_t<cluster_tag_t>& cluster_assignments,
                    sequence_data_t<short> tfbs_start_positions) {
                std::list<const index_t*> binding_sites;
                // find all binding sites
                for (Indexer::sampling_iterator it = indexer.sampling_begin();
                     it != indexer.sampling_end(); it++) {
                        if (tfbs_start_positions[**it] == 1) {
                                binding_sites.push_back(*it);
                        }
                }
                // iterate over binding sites
                for (std::list<const index_t*>::const_iterator it = binding_sites.begin();
                     it != binding_sites.end(); it++) {
                        // if there still is a binding site
                        if (tfbs_start_positions[**it] == 1) {
                                cluster_tag_t tag = cluster_assignments[**it];
                                tfbs_start_positions[**it] = 0;
                                // find sites with the same cluster assignment
                                for (std::list<const index_t*>::const_iterator is = it;
                                     is != binding_sites.end(); is++) {
                                        if (tfbs_start_positions[**is] == 1 && cluster_assignments[**is] == tag) {
                                                tfbs_start_positions[**is] = 0;
                                                _edges[edge_t(**it, **is)]++;
                                        }
                                }
                        }
                }
        }
        void print() const {
                for (map_t::const_iterator it = _edges.begin();
                     it != _edges.end(); it++) {
                        if ((*it).second >= 10) {
                                std::cout << (*it).first.index1 << ":" << (*it).first.index2 << " -> " << (*it).second << std::endl;
                        }
                }
        }

protected:
        map_t _edges;
};

#endif /* DPM_TFBS_GRAPH_HH */
