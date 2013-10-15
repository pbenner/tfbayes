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
#include <ostream>
#include <map>

#include <boost/format.hpp>
#include <boost/unordered_map.hpp>

#include <tfbayes/fasta/fasta.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/uipac/code.hh>
#include <tfbayes/uipac/nucleotide-sequence.hh>
#include <tfbayes/utility/linalg.hh>
#include <tfbayes/utility/strtools.hh>
#include <tfbayes/alignment/sequence.hh>
#include <tfbayes/dpm/index.hh>

typedef seq_index_t alignment_index_t;

template <typename CODE_TYPE = alphabet_code_t>
class alignment_t : std::matrix<CODE_TYPE> {
public:
        using std::matrix<CODE_TYPE>::begin;
        using std::matrix<CODE_TYPE>::end;
        // Typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_map<std::string, pt_node_t::id_t> taxon_map_t;
        typedef typename std::matrix<CODE_TYPE>::iterator iterator;
        typedef typename std::matrix<CODE_TYPE>::const_iterator const_iterator;

        // Constructors
        ////////////////////////////////////////////////////////////////////////
        alignment_t(size_t length,
                    const pt_root_t& tree,
                    CODE_TYPE init = -1,
                    alphabet_t alphabet = nucleotide_alphabet_t())
                : std::matrix<CODE_TYPE>(length, tree.n_leaves, init),
                  _n_species (tree.n_leaves),
                  _length    (length),
                  _alphabet  (alphabet) {
                // initialize taxon map from leaf names
                init_taxon_map(tree);
        }
        alignment_t(const std::matrix<CODE_TYPE> sequences,
                    const taxon_map_t& taxon_map,
                    alphabet_t alphabet = nucleotide_alphabet_t())
                : std::matrix<CODE_TYPE>(sequences),
                  _n_species (taxon_map.size()),
                  // the length is initialized later
                  _length    (0),
                  _taxon_map (taxon_map),
                  _alphabet  (alphabet) {
                // check that all lengths are consistent
                init_alignment(sequences);
        }
        alignment_t(const std::string& filename,
                    boost::optional<const pt_root_t&> tree = boost::optional<const pt_root_t&>(),
                    alphabet_t alphabet = nucleotide_alphabet_t())
                : std::matrix<CODE_TYPE>(),
                  // we have as many sequences in this alignment as
                  // there are leaves in the tree
                  _n_species (0),
                  // the length is initialized later
                  _length    (0),
                  _alphabet  (alphabet)
                {
                // fill taxon map seperately since some species might
                // not be present in the alignment
                tree ? init_taxon_map(*tree   ):
                       init_taxon_map(filename);
                // with the taxon map initialized we know the number
                // of species in the alignment
                _n_species = _taxon_map.size();
                // parse fasta file
                std::matrix<CODE_TYPE> tmp = parse_fasta(filename);
                // check that all lengths are consistent
                init_alignment(tmp);
        }
        alignment_t(const alignment_t& alignment)
                : std::matrix<CODE_TYPE>(alignment),
                  _n_species (alignment._n_species),
                  _length    (alignment._length),
                  _taxon_map (alignment._taxon_map),
                  _alphabet  (alignment._alphabet) {
        }
        virtual ~alignment_t() {
        }

        friend void swap(alignment_t& first, alignment_t&second) {
                std::swap(static_cast<std::matrix<CODE_TYPE>&>(first),
                          static_cast<std::matrix<CODE_TYPE>&>(second));
                std::swap(first._length,    second._length);
                std::swap(first._n_species, second._n_species);
                std::swap(first._taxon_map, second._taxon_map);
                std::swap(first._alphabet,  second._alphabet);
        }

        // Operators
        ////////////////////////////////////////////////////////////////////////
        alignment_t& operator=(const alignment_t& alignment) {
                alignment_t tmp(alignment);
                swap(*this, tmp);
                return *this;
        }
        // use the access operator from base class to read/write columns
        using std::matrix<CODE_TYPE>::operator[];
        // define new access operators to access individual cells
        const CODE_TYPE& operator[](const alignment_index_t& index) const {
                return std::matrix<CODE_TYPE>::operator[](index[1])[index[0]];
        }
              CODE_TYPE& operator[](const alignment_index_t& index) {
                return std::matrix<CODE_TYPE>::operator[](index[1])[index[0]];
        }
        // and to obtain full sequences for one species
        sequence_t<CODE_TYPE> operator[](const std::string& taxon) const {
                sequence_t<CODE_TYPE> sequence(length(), alphabet());
                pt_node_t::id_t id = taxon_map(taxon);
                for (const_iterator it = begin(); it != end() && id != -1; it++) {
                        sequence[it-begin()] = (*it)[id];
                }
                return sequence;
        }

        // Access Methods
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
        const alphabet_t& alphabet() const {
                return _alphabet;
        }
        // Methods
        ////////////////////////////////////////////////////////////////////////
        template<size_t ALPHABET_SIZE>
        std::vector<double> scan(const pt_root_t& tree, std::matrix<double>& counts) {
                std::vector<double> result(length(), 0);
                vector<exponent_t<CODE_TYPE, ALPHABET_SIZE> > exponents;

                for (size_t j = 0; j < counts.size(); j++) {
                        exponent_t<CODE_TYPE, ALPHABET_SIZE> tmp
                                (counts[j].begin(), counts[j].end());
                        exponents.push_back(tmp);
                }

                for (iterator it = begin(); it != end(); it++) {
                        result[it - begin()] = 0;
                        if (it - begin() + counts.size() > length()) {
                                // do not exit the loop here so that every
                                // position is initialized
                                continue;
                        }
                        for (iterator is(it); is < it + counts.size(); is++) {
                                result[it - begin()] += pt_marginal_likelihood<CODE_TYPE, ALPHABET_SIZE>(
                                        tree, *is, exponents[is-it]);
                        }
                }
                return result;
        }
        template<size_t ALPHABET_SIZE>
        std::vector<double> marginal_likelihood(const pt_root_t& tree, const std::vector<double>& prior) {
                exponent_t<CODE_TYPE, ALPHABET_SIZE> alpha(prior.begin(), prior.end());
                vector<double> result;

                /* go through the alignment and compute the marginal
                 * likelihood for each position */
                for (iterator it = begin(); it != end(); it++) {
                        result.push_back(pt_marginal_likelihood<CODE_TYPE, ALPHABET_SIZE>
                                         (tree, *it, alpha));
                }
                return result;
        }
protected:
        // Methods for Initializations
        ////////////////////////////////////////////////////////////////////////
        void init_taxon_map(const pt_root_t& tree) {
                for (pt_node_t::leaves_t::const_iterator it = tree.leaves.begin();
                     it != tree.leaves.end(); it++) {
                        _taxon_map[(*it)->name] = (*it)->id;
                }
        }
        void init_taxon_map(const std::string& filename) {
                /* this automatically parses the fasta file format */
                FastaParser parser(filename);
                std::set<std::string> species;

                while (parser) {
                        if (parser() == "")
                                continue;
                        if (parser.taxon() == "")
                                continue;
                        species.insert(parser.taxon());
                }
                std::vector<std::string> tmp(species.begin(), species.end());
                for (size_t i = 0; i < tmp.size(); i++) {
                        _taxon_map[tmp[i]] = i;
                }
        }
        std::matrix<CODE_TYPE> parse_fasta(const std::string& filename) {
                FastaParser parser(filename);

                std::matrix<CODE_TYPE> tmp(n_species(), 0);
                std::string line;

                while (parser) {
                        line = parser();
                        if (line == "")
                                continue;
                        if (parser.description().size() == 0) {
                                std::cerr << "Warning: sequence without description found... skipping sequence."
                                          << std::endl;
                                continue;
                        }
                        pt_node_t::id_t id = taxon_map(parser.taxon());
                        if (id != -1) {
                                sequence_t<CODE_TYPE> sequence(line, alphabet());
                                // multiple entries should be ignored!
                                if (tmp[id].size() > 0) {
                                        std::cerr << boost::format("Warning: taxon `%s' appeared more than once in the alignment... "
                                                                   "skipping sequence.") % parser.taxon()
                                                  << std::endl;
                                }
                                else {
                                        tmp[id].insert(tmp[id].end(), sequence.begin(), sequence.end());
                                }
                        }
                        else {
                                std::cerr << boost::format("Warning: taxon `%s' not found in the phylogenetic tree.") % parser.taxon()
                                          << std::endl;
                        }
                }
                return tmp;
        }
        void init_alignment(const std::matrix<CODE_TYPE>& sequences) {
                std::matrix<CODE_TYPE> tmp(sequences);
                // initialize length
                for (typename std::matrix<CODE_TYPE>::const_iterator it = tmp.begin();
                     it != tmp.end(); it++) {
                        // update length if necessary
                        if (it->size() > length()) {
                                _length = it->size();
                        }
                }
                // initialize data
                for (taxon_map_t::const_iterator it = taxon_map().begin();
                     it != taxon_map().end(); it++) {
                        if (tmp[it->second].size() == 0) {
                                std::cerr << boost::format("Warning: taxon `%s' has sequence length zero") % it->first
                                          << std::endl;
                                // fill with missing data
                                tmp[it->second] = sequence_t<CODE_TYPE>(length(), alphabet());
                        }
                        else
                        if (tmp[it->second].size() != length()) {
                                std::cerr << boost::format("Error: taxon `%s' has inconsistent sequence length: ") % it->first
                                          << boost::format("%d instead of %d") % tmp[it->second].size() % length()
                                          << std::endl;
                                exit(EXIT_FAILURE);
                        }
                }
                // assign contents of tmp to this object
                std::matrix<CODE_TYPE>::operator=(tmp.transpose());
        }
        // Fields
        ////////////////////////////////////////////////////////////////////////
        // number of species
        size_t _n_species;
        // length of the nucleotide sequence
        size_t _length;
        // map taxon name to nodes in the phylogenetic tree
        taxon_map_t _taxon_map;
        // alphabet (i.e. nucleotides or amino acids)
        alphabet_t _alphabet;
};


template <typename CODE_TYPE = alphabet_code_t>
class alignment_set_t : public std::vector<alignment_t<CODE_TYPE> > {
public:
        using std::vector<alignment_t<CODE_TYPE> >::push_back;

        // Typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_map<std::string, pt_node_t::id_t> taxon_map_t;

        // Constructors
        ////////////////////////////////////////////////////////////////////////
        alignment_set_t() { };
        alignment_set_t(const std::string& filename,
                        boost::optional<const pt_root_t&> tree = boost::optional<const pt_root_t&>(),
                        alphabet_t alphabet = nucleotide_alphabet_t()) {
                /* first check what species are available */
                taxon_map_t taxon_map = tree ?
                        create_taxon_map(*tree   ):
                        create_taxon_map(filename);

                /* number of species */
                size_t n = taxon_map.size();

                /* exit if this file is empty */
                if (n == 0) return;

                /* this automatically parses the fasta file format */
                FastaParser parser(filename);

                /* storage for a single line */
                std::string line;

                /* remember which species already occured */
                std::set<std::string> occurred;

                /* current alignment */
                std::matrix<CODE_TYPE> sequences(n, 0);

                /* the fasta parser returns a single line for each
                 * sequence */
                while (parser) {
                        line = parser();
                        if (line == "")
                                continue;
                        if (parser.description().size() == 0) {
                                std::cerr << "Warning: sequence without description found... skipping."
                                          << std::endl;
                                continue;
                        }
                        if (occurred.find(parser.taxon()) != occurred.end()) {
                                // push alignment
                                push_back(alignment_t<CODE_TYPE>(sequences, taxon_map, alphabet));
                                // reset occurrences
                                occurred.clear();
                                // start new alignment
                                sequences = std::matrix<CODE_TYPE>(n, 0);
                        }
                        occurred.insert(parser.taxon());
                        sequences[taxon_map[parser.taxon()]] = nucleotide_sequence_t<CODE_TYPE>(line);
                }
                push_back(alignment_t<CODE_TYPE>(sequences, taxon_map, alphabet));
        }

protected:
        // Methods for Initializations
        ////////////////////////////////////////////////////////////////////////
        taxon_map_t create_taxon_map(const pt_root_t& tree) {
                taxon_map_t taxon_map;

                for (pt_node_t::leaves_t::const_iterator it = tree.leaves.begin();
                     it != tree.leaves.end(); it++) {
                        taxon_map[(*it)->name] = (*it)->id;
                }
                return taxon_map;
        }
        taxon_map_t create_taxon_map(const std::string& filename) {
                /* this automatically parses the fasta file format */
                FastaParser parser(filename);
                std::set<std::string> species;
                taxon_map_t taxon_map;

                while (parser) {
                        if (parser() == "")
                                continue;
                        if (parser.taxon() == "")
                                continue;
                        species.insert(parser.taxon());
                }
                std::vector<std::string> tmp(species.begin(), species.end());
                for (size_t i = 0; i < tmp.size(); i++) {
                        taxon_map[tmp[i]] = i;
                }
                return taxon_map;
        }
};

template <typename CODE_TYPE>
class print_alignment_t {
public:
        print_alignment_t(const alignment_t<CODE_TYPE>& alignment)
                : _alignment(alignment)
                { }

        virtual const alignment_t<CODE_TYPE>& operator()() const {
                return _alignment;
        }
protected:
        const alignment_t<CODE_TYPE>& _alignment;
};

template <typename CODE_TYPE>
class print_alignment_pretty_t : public print_alignment_t<CODE_TYPE> {
public:
        print_alignment_pretty_t(const alignment_t<CODE_TYPE>& alignment)
                : print_alignment_t<CODE_TYPE>(alignment)
                { }
};
template <typename CODE_TYPE>
print_alignment_pretty_t<CODE_TYPE> print_alignment_pretty(alignment_t<CODE_TYPE> alignment)
{
        return print_alignment_pretty_t<CODE_TYPE>(alignment);
}

template <typename CODE_TYPE>
class print_alignment_fasta_t : public print_alignment_t<CODE_TYPE> {
public:
        print_alignment_fasta_t(const alignment_t<CODE_TYPE>& alignment)
                : print_alignment_t<CODE_TYPE>(alignment)
                { }
};
template <typename CODE_TYPE>
print_alignment_fasta_t<CODE_TYPE> print_alignment_fasta(alignment_t<CODE_TYPE> alignment)
{
        return print_alignment_fasta_t<CODE_TYPE>(alignment);
}

template <typename CODE_TYPE>
std::ostream& operator<< (std::ostream& o, const print_alignment_pretty_t<CODE_TYPE>& alignment_container)
{
        const alignment_t<CODE_TYPE>& alignment(alignment_container());

        for (typename alignment_t<CODE_TYPE>::taxon_map_t::const_iterator it = alignment.taxon_map().begin();
             it != alignment.taxon_map().end(); it++) {

                o << boost::format("%10s (%2d): ") % it->first % it->second
                  << alignment[it->first]
                  << std::endl;
        }
        return o;
}

template <typename CODE_TYPE>
std::ostream& operator<< (std::ostream& o, const print_alignment_fasta_t<CODE_TYPE>&  alignment_container)
{
        const alignment_t<CODE_TYPE>& alignment(alignment_container());

        for (typename alignment_t<CODE_TYPE>::taxon_map_t::const_iterator it = alignment.taxon_map().begin();
             it != alignment.taxon_map().end(); it++) {
                std::stringstream ss;
                ss << alignment[it->first];
                o << boost::format(">%s\n") % it->first
                  << split_string(ss.str(), 80)
                  << std::endl;
        }
        return o;
}

#endif /* ALIGNMENT_HH */
