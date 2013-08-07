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
#include <vector>

#include <boost/unordered_map.hpp>

#include <tfbayes/fasta/fasta.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/uipac/code.hh>
#include <tfbayes/uipac/nucleotide-sequence.hh>
#include <tfbayes/utility/strtools.hh>

template <typename CODE_TYPE>
class alignment_t : std::vector<nucleotide_sequence_t<CODE_TYPE> > {
public:
        // Typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef std::vector<nucleotide_sequence_t<CODE_TYPE> > alignment_ancestor_t;
        typedef boost::unordered_map<std::string, pt_node_t::id_t> taxon_map_t;

        // Constructors
        ////////////////////////////////////////////////////////////////////////
        alignment_t(const alignment_ancestor_t sequences, const taxon_map_t& taxon_map)
                : alignment_ancestor_t(sequences),
                  _n_species(sequences.size()),
                  _taxon_map(taxon_map) {
                check_lengths();
        }
        alignment_t(const size_t length, CODE_TYPE init, const pt_root_t* tree)
                : alignment_ancestor_t(tree->n_leaves, nucleotide_sequence_t<CODE_TYPE>()),
                  _n_species(tree->n_leaves), _length(length) {
                for (pt_node_t::leaves_t::const_iterator it = tree->leaves.begin(); it != tree->leaves.end(); it++) {
                        pt_leaf_t* leaf = *it;
                        _taxon_map[leaf->name]     = leaf->id;
                        this->operator[](leaf->id) = nucleotide_sequence_t<CODE_TYPE>(length, init);
                }
        }
        alignment_t(const char* filename, const pt_root_t* tree)
                : alignment_ancestor_t(tree->n_leaves, nucleotide_sequence_t<CODE_TYPE>()) {

                FastaParser parser(filename);
                std::string sequence;

                // we have as many sequences in this alignment as
                // there are leaves in the tree
                _n_species = tree->n_leaves;
                // fill taxon map seperately since some species might
                // not be present in the alignment
                for (pt_node_t::leaves_t::const_iterator it = tree->leaves.begin();
                     it != tree->leaves.end(); it++) {
                        _taxon_map[(*it)->name] = (*it)->id;
                }
                // parse fasta file
                while ((sequence = parser.read_sequence()) != "") {
                        std::string taxon  = token(parser.description()[0], '.')[0];
                        pt_node_t::id_t id = tree->get_leaf_id(taxon);
                        if (id != -1) {
                                // there might be multiple entries for
                                // each species, join all entries into
                                // a single sequence
                                operator[](id) = nucleotide_sequence_t<CODE_TYPE>(operator[](id), sequence);
                        }
                        else {
                                std::cerr << "Warning: taxon `"
                                          << taxon
                                          << "' not found in the phylogenetic tree."
                                          << std::endl;
                        }
                }
                check_lengths();
        }
        alignment_t(const alignment_t& alignment)
                : alignment_ancestor_t() {
                operator=(alignment);
        }
        virtual ~alignment_t() {
        }

        // Operators
        ////////////////////////////////////////////////////////////////////////
        alignment_t& operator=(const alignment_t& alignment) {
                alignment_ancestor_t::operator=(alignment);
                _length    = alignment._length;
                _n_species = alignment._n_species;
                _taxon_map = alignment._taxon_map;
                return *this;
        }
        using alignment_ancestor_t::operator[];
        const nucleotide_sequence_t<CODE_TYPE>& operator[](const std::string& taxon) const {
                return operator[](taxon_map(taxon));
        }
        nucleotide_sequence_t<CODE_TYPE>& operator[](const std::string& taxon) {
                return operator[](taxon_map(taxon));
        }

        // Iterator
        ////////////////////////////////////////////////////////////////////////
        class iterator
        {
        public:
                iterator(size_t position,
                         size_t length,
                         size_t n_species,
                         const alignment_ancestor_t& alignment,
                         const taxon_map_t& taxon_map)
                        : _position(position),
                          _length(length),
                          _n_species(n_species),
                          _alignment(alignment),
                          _observations(n_species, 0) {
                        fill_observations();
                }
                iterator(const iterator& it)
                        : _position(0),
                          _length(it._length),
                          _n_species(it._n_species),
                          _alignment(it._alignment),
                          _observations(it._observations)
                {
                        operator=(it);
                }
                iterator& operator=(const iterator& it)
                {
                        _position     = it._position;
                        _observations = it._observations;

                        return *this;
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
                                fill_observations();
                        }
                        return *this;
                }
                iterator& operator--(int _) {
                        if (_position > 0) {
                                _position--;
                                fill_observations();
                        }
                        return *this;
                }
                const std::vector<CODE_TYPE>& operator*() const {
                        return _observations;
                }
                const std::vector<CODE_TYPE>* operator->() const {
                        return &operator*();
                }
        protected:
                void fill_observations() {
                        if (position() < _length) {
                                for (size_t i = 0; i < _n_species; i++) {
                                        if (_alignment[i].size() > 0) {
                                                _observations[i] = _alignment[i][position()];
                                        }
                                        else {
                                                _observations[i] = code_nucleotide<CODE_TYPE>('-');
                                        }
                                }
                        }
                }

                size_t _position;
                const size_t _length;
                const size_t _n_species;
                const alignment_ancestor_t& _alignment;
                std::vector<CODE_TYPE> _observations;
        };
        virtual iterator begin() const {
                return iterator(0, _length, _n_species, *this, _taxon_map);
        }
        virtual iterator end() const {
                return iterator(_length, _length, _n_species, *this, _taxon_map);
        }
        // Methods
        ////////////////////////////////////////////////////////////////////////
        const size_t& n_species() const {
                return _n_species;
        }
        const size_t& length() const {
                return _length;
        }
        const taxon_map_t& taxon_map() const {
                return _taxon_map;
        }
        pt_node_t::id_t taxon_map(const std::string& taxon) const {
                taxon_map_t::const_iterator it = _taxon_map.find(taxon);
                if (it != _taxon_map.end()) {
                        return it->second;
                }
                return -1;
        }

protected:
        // number of species
        size_t _n_species;
        // length of the nucleotide sequence
        size_t _length;
        // map taxon name to nodes in the phylogenetic tree
        taxon_map_t _taxon_map;

        // private methods
        void check_lengths() {
                // find a sequence that has length greater than zero
                // to set the length of this alignment
                _length = 0;
                for (typename alignment_ancestor_t::const_iterator it = alignment_ancestor_t::begin();
                     it != alignment_ancestor_t::end(); it++) {
                        if (it->size() > 0) {
                                _length = it->size();
                                break;
                        }
                }
                // check that every sequence has either length zero or
                // length equal to the length of this alignment
                for (taxon_map_t::const_iterator it = _taxon_map.begin();
                     it != _taxon_map.end(); it++) {
                        if (operator[](it->second).size() == 0) {
                                std::cerr << "Warning: nucleotide sequence for taxon `"
                                          << it->first
                                          << "' has length zero."
                                          << std::endl;
                        }
                        if(operator[](it->second).size() != 0 &&
                           operator[](it->second).size() != _length) {
                                std::cerr << "Warning: nucleotide sequence for taxon `"
                                          << it->first
                                          << "' has invalid length."
                                          << std::endl;
                                exit(EXIT_FAILURE);
                        }
                }
        }
};


template <typename CODE_TYPE>
class alignment_set_t : public std::vector<alignment_t<CODE_TYPE> > {
public:
        // Typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_map<std::string, pt_node_t::id_t> taxon_map_t;

        // Constructors
        ////////////////////////////////////////////////////////////////////////
        alignment_set_t() { };
        alignment_set_t(const char* filename) {
                read_fasta(filename);
        }

private:
        taxon_map_t get_species_from_fasta(const char* filename) {
                /* this automatically parses the fasta file format */
                FastaParser parser(filename);
                std::set<std::string> species;
                taxon_map_t taxon_map;

                while (parser.read_sequence() != "") {
                        assert(parser.description().size()    > 0);
                        assert(parser.description()[0].size() > 0);
                        species.insert(parser.description()[0]);
                }
                std::vector<std::string> tmp(species.begin(), species.end());
                for (size_t i = 0; i < tmp.size(); i++) {
                        taxon_map[tmp[i]] = i;
                }
                return taxon_map;
        }

        void read_fasta(const char* filename) {
                /* first check what species are available */
                taxon_map_t taxon_map = get_species_from_fasta(filename);

                /* number of species */
                size_t n = taxon_map.size();

                /* exit if this file is empty */
                if (n == 0) return;

                /* this automatically parses the fasta file format */
                FastaParser parser(filename);

                /* storage for a single line */
                std::string line;

                /* remember which species already occured */
                std::map<std::string, bool> occurred;

                /* current alignment */
                std::vector<nucleotide_sequence_t<CODE_TYPE> > sequences(n, nucleotide_sequence_t<CODE_TYPE>());

                /* the fasta parser returns a single line for each
                 * sequence */
                while ((line = parser.read_sequence()) != "") {
                        const std::string& name = parser.description()[0];
                        if (occurred[name]) {
                                // push alignment
                                this->push_back(alignment_t<CODE_TYPE>(sequences, taxon_map));
                                // reset occurrences
                                occurred = std::map<std::string, bool>();
                                // start new alignment
                                sequences = std::vector<nucleotide_sequence_t<CODE_TYPE> >(n, nucleotide_sequence_t<CODE_TYPE>());
                        }
                        occurred[name] = true;
                        sequences[taxon_map[name]] = nucleotide_sequence_t<CODE_TYPE>(line);
                }
                this->push_back(alignment_t<CODE_TYPE>(sequences, taxon_map));
        }
};

#endif /* ALIGNMENT_HH */
