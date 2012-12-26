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

#ifndef ALIGNMENT_HH
#define ALIGNMENT_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <set>
#include <string>

#include <boost/unordered_map.hpp>

#include <tfbayes/fasta/fasta.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/uipac/nucleotide-sequence.hh>
#include <tfbayes/utility/strtools.hh>

template <typename CODE_TYPE>
class alignment_t : boost::unordered_map<std::string, nucleotide_sequence_t<CODE_TYPE> > {
public:
        // Typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_map<std::string, nucleotide_sequence_t<CODE_TYPE> > alignment_ancestor_t;
        typedef boost::unordered_map<std::string, pt_node_t::id_t> taxon_map_t;

        // Constructors
        ////////////////////////////////////////////////////////////////////////
        alignment_t(const size_t length, pt_root_t* tree)
                : alignment_ancestor_t(),
                  _length(length),
                  _tree(tree->clone()) {

                pt_node_t::nodes_t nodes = _tree->get_nodes();

                for (pt_node_t::nodes_t::const_iterator it = nodes.begin(); it != nodes.end(); it++) {
                        pt_node_t* node = *it;
                        _taxon_map[node->name] = node->id;
                        this->operator[](node->name) = nucleotide_sequence_t<CODE_TYPE>(length);
                }
                _size = nodes.size();
        }
        alignment_t(const char* filename, pt_root_t* tree)
                : alignment_ancestor_t(),
                  _tree(tree->clone()) {

                FastaParser parser(filename);
                std::string sequence;

                while ((sequence = parser.read_sequence()) != "") {
                        std::string taxon  = token(parser.description()[0], '.')[0];
                        pt_node_t::id_t id = _tree->get_id(taxon);
                        if (id != -1) {
                                _taxon_map[taxon] = id;
                                operator[](taxon) = nucleotide_sequence_t<CODE_TYPE>(sequence);
                        }
                        else {
                                std::cerr << "Warning: taxon `"
                                          << taxon
                                          << "' not found in the phylogenetic tree."
                                          << std::endl;
                        }
                }
                if (alignment_ancestor_t::size() > 0) {
                        _length = alignment_ancestor_t::begin()->second.size();
                }
                else {
                        _length = 0;
                }
                _size = _taxon_map.size();
        }
        alignment_t(const alignment_t& alignment)
                : alignment_ancestor_t(),
                  _tree(NULL) {
                operator=(alignment);
        }
        virtual ~alignment_t() {
                if (_tree) {
                        _tree->destroy();
                }
        }
        // Operators
        ////////////////////////////////////////////////////////////////////////
        alignment_t& operator=(const alignment_t& alignment) {
                alignment_ancestor_t::operator=(alignment);
                if (_tree) {
                        _tree->destroy();
                }
                _tree      = alignment._tree->clone();
                _taxon_map = alignment._taxon_map;
                return *this;
        }
        using alignment_ancestor_t::operator[];

        // Iterator
        ////////////////////////////////////////////////////////////////////////
        class iterator
        {
        public:
                iterator(size_t position,
                         size_t length,
                         pt_root_t *tree,
                         const alignment_ancestor_t& alignment,
                         const taxon_map_t& taxon_map)
                        : _position(position),
                          _tree(tree->clone()),
                          _length(length),
                          _alignment(alignment),
                          _taxon_map(taxon_map)
                { }
                iterator(const iterator& it)
                        : _position(0),
                          _tree(NULL),
                          _length(it._length),
                          _alignment(it._alignment),
                          _taxon_map(it._taxon_map)
                {
                        operator=(it);
                }
                ~iterator()
                {
                        if (_tree) {
                                _tree->destroy();
                        }
                }
                iterator& operator=(const iterator& it)
                {
                        _position = it._position;
                        if (_tree) {
                                _tree->destroy();
                        }
                        _tree = it._tree->clone();

                        return *this;
                }
                void apply(pt_root_t* tree) {
                }
                const size_t& position() const {
                        return _position;
                }
                bool operator==(const iterator& it) const {
                        return position() == it.position();
                }
                bool operator!=(const iterator& it) const {
                        return position() != it.position();
                }
                iterator& operator++(int _) {
                        if (_position < _length) {
                                _position++;
                        }
                        return *this;
                }
                iterator& operator--(int _) {
                        if (_position > 0) {
                                _position--;
                        }
                        return *this;
                }
                const pt_root_t& operator*() {
                        for (taxon_map_t::const_iterator it = _taxon_map.begin();
                             it != _taxon_map.end(); it++) {
                                const std::string taxon  = it->first;
                                const pt_node_t::id_t id = it->second;
                                _tree->node_map[id]->x = _alignment.find(taxon)->second[position()];
                        }
                        return *_tree;
                }
                const pt_root_t* operator->() {
                        return &operator*();
                }
        protected:
                size_t _position;
                pt_root_t* _tree;
                const size_t _length;
                const alignment_ancestor_t& _alignment;
                const taxon_map_t& _taxon_map;
        };
        virtual iterator begin() const {
                return iterator(0, _length, _tree, *this, _taxon_map);
        }
        virtual iterator end() const {
                return iterator(_length, _length, _tree, *this, _taxon_map);
        }
        // Methods
        ////////////////////////////////////////////////////////////////////////
        const size_t& size() const {
                return _size;
        }
        const size_t& length() const {
                return _length;
        }

protected:
        // number of species
        size_t _size;
        // length of the nucleotide sequence
        size_t _length;
        // map taxon name to nodes in the phylogenetic tree
        taxon_map_t _taxon_map;
        // phylogenetic tree
        pt_root_t* _tree;
};

#endif /* ALIGNMENT_HH */
