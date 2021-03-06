/* Copyright (C) 2012-2013 Philipp Benner
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

#ifndef __TFBAYES_ALIGNMENT_ALIGNMENT_HH__
#define __TFBAYES_ALIGNMENT_ALIGNMENT_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>    // std::random_shuffle
#include <set>
#include <string>
#include <vector>
#include <ostream>
#include <map>

#include <boost/format.hpp>
#include <boost/unordered_map.hpp>

#include <tfbayes/fasta/fasta.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/linalg.hh>
#include <tfbayes/utility/strtools.hh>
#include <tfbayes/alignment/sequence.hh>
#include <tfbayes/dpm/index.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

template <typename AC = alphabet_code_t>
class alignment_t : std::matrix<AC> {
public:
        using std::matrix<AC>::begin;
        using std::matrix<AC>::end;
        // Typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef std::matrix<AC> base_t;
        typedef boost::unordered_map<std::string, pt_node_t::id_t> taxon_map_t;
        typedef typename base_t::iterator iterator;
        typedef typename base_t::const_iterator const_iterator;

        // Constructors
        ////////////////////////////////////////////////////////////////////////
        alignment_t(size_t length,
                    const pt_root_t& tree,
                    AC init = -1,
                    alphabet_t alphabet = nucleotide_alphabet_t())
                : base_t     (length, tree.n_leaves, init),
                  _n_species (tree.n_leaves),
                  _length    (length),
                  _alphabet  (alphabet) {
                // initialize taxon map from leaf names
                init_taxon_map(tree);
        }
        alignment_t(const std::matrix<AC>& sequences,
                    const taxon_map_t& taxon_map,
                    alphabet_t alphabet = nucleotide_alphabet_t(),
                    bool verbose = false)
                : base_t     (sequences),
                  _n_species (taxon_map.size()),
                  // the length is initialized later
                  _length    (0),
                  _taxon_map (taxon_map),
                  _alphabet  (alphabet) {
                // check that all lengths are consistent
                init_alignment(sequences, verbose);
        }
        alignment_t(const std::matrix<AC>& sequences,
                    const pt_root_t& tree,
                    alphabet_t alphabet = nucleotide_alphabet_t(),
                    bool verbose = false)
                : base_t     (sequences),
                  _n_species (tree.n_leaves),
                  _length    (0),
                  _alphabet  (alphabet) {
                init_taxon_map(tree);
                // check that all lengths are consistent
                init_alignment(sequences, verbose);
        }
        alignment_t(const std::string& filename,
                    boost::optional<const pt_root_t&> tree = boost::optional<const pt_root_t&>(),
                    alphabet_t alphabet = nucleotide_alphabet_t(),
                    bool verbose = false)
                : base_t(),
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
                std::matrix<AC> tmp = parse_fasta(filename, verbose);
                // check that all lengths are consistent
                init_alignment(tmp, verbose);
        }
        alignment_t(const alignment_t& alignment)
                : base_t(alignment),
                  _n_species (alignment._n_species),
                  _length    (alignment._length),
                  _taxon_map (alignment._taxon_map),
                  _alphabet  (alignment._alphabet) {
        }
        virtual ~alignment_t() {
        }

        friend void swap(alignment_t& first, alignment_t&second) {
                std::swap(static_cast<base_t&>(first),
                          static_cast<base_t&>(second));
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
        using base_t::operator[];
        // define new access operators to access individual cells
        const AC& operator[](const index_t& index) const {
                return base_t::operator[](index[1])[index[0]];
        }
              AC& operator[](const index_t& index) {
                return base_t::operator[](index[1])[index[0]];
        }
        alignment_t<AC> operator[](const range_t& range) const {
                assert(typeid(range.index()) == typeid(index_t));
                base_t sequences;
                const index_t& index = range.index();
                for (size_t i = 0; i < range.length(); i++) {
                        sequences.push_back(base_t::operator[](index[0]+i));
                }
                return alignment_t<AC>(sequences.transpose(), taxon_map(), alphabet());
        }
        // and to obtain full sequences for one species
        sequence_t<AC> operator[](const std::string& taxon) const {
                sequence_t<AC> sequence(length(), alphabet());
                pt_node_t::id_t id = taxon_map(taxon);
                for (const_iterator it = begin(); it != end() && id != -1; it++) {
                        sequence[it-begin()] = (*it)[id];
                }
                return sequence;
        }
        void operator+=(const alignment_t& alignment) {
                assert(taxon_map() == alignment.taxon_map());
                assert(alphabet()  == alignment.alphabet());
                for (const_iterator it = alignment.begin();
                     it != alignment.end(); it++) {
                        base_t::push_back(*it);
                }
                _length += alignment.length();
        }
        void shuffle() {
                // this is not thread safe!
                std::random_shuffle(begin(), end());
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
protected:
        // Methods for Initializations
        ////////////////////////////////////////////////////////////////////////
        void init_taxon_map(const pt_root_t& tree) {
                for (pt_node_t::leaves_t::const_iterator it = tree.begin_leaves();
                     it != tree.end_leaves(); it++) {
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
        std::matrix<AC> parse_fasta(const std::string& filename, bool verbose = false) {
                FastaParser parser(filename);

                std::matrix<AC> tmp(n_species(), 0);
                std::string line;

                for (size_t i = 1; parser; i++) {
                        if (verbose) {
                                std::cerr << boost::format("Reading sequence %d...") % i
                                          << std::endl;
                        }
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
                                sequence_t<AC> sequence(line, alphabet());
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
        void init_alignment(const std::matrix<AC>& sequences, bool verbose = false) {
                std::matrix<AC> tmp(sequences);
                // initialize length
                for (typename base_t::const_iterator it = tmp.begin();
                     it != tmp.end(); it++) {
                        // update length if necessary
                        if (it->size() > length()) {
                                _length = it->size();
                        }
                }
                if (verbose) {
                        std::cerr << "-> alignment length is " << _length << std::endl;
                }
                // initialize data
                for (taxon_map_t::const_iterator it = taxon_map().begin();
                     it != taxon_map().end(); it++) {
                        if (tmp[it->second].size() == 0) {
                                std::cerr << boost::format("Warning: taxon `%s' has sequence length zero (regarding as missing data)") % it->first
                                          << std::endl;
                                // fill with missing data
                                tmp[it->second] = sequence_t<AC>(length(), alphabet());
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
                base_t::operator=(tmp.transpose());
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


template <typename AC = alphabet_code_t>
class alignment_set_t : public std::vector<alignment_t<AC> > {
public:
        using std::vector<alignment_t<AC> >::push_back;

        // Typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef std::vector<alignment_t<AC> > base_t;
        typedef boost::unordered_map<std::string, pt_node_t::id_t> taxon_map_t;

        // Constructors
        ////////////////////////////////////////////////////////////////////////
        alignment_set_t()
                : base_t()
                { };
        alignment_set_t(const std::string& filename,
                        boost::optional<const pt_root_t&> tree = boost::optional<const pt_root_t&>(),
                        alphabet_t alphabet = nucleotide_alphabet_t(),
                        bool verbose = false)
                : base_t() {
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
                std::matrix<AC> sequences(n, 0);

                /* the fasta parser returns a single line for each
                 * sequence */
                for (size_t i = 1; parser;) {
                        line = parser();
                        if (parser.description().size() == 0) {
                                std::cerr << "Warning: sequence without description found... skipping."
                                          << std::endl;
                                continue;
                        }
                        if (occurred.find(parser.taxon()) != occurred.end() ||
                            parser.description()[0] == ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>") {
                                // push alignment
                                push_back(alignment_t<AC>(sequences, taxon_map, alphabet, verbose));
                                if (verbose) {
                                        std::cerr << boost::format("Finished parsing alignment %d...") % i++
                                                  << std::endl;
                                }
                                // reset occurrences
                                occurred.clear();
                                // start new alignment
                                sequences = std::matrix<AC>(n, 0);
                                continue;
                        }
                        if (line == "") {
                                std::cerr << "Warning: empty sequence found... skipping."
                                          << std::endl;
                                continue;
                        }
                        taxon_map_t::const_iterator it = taxon_map.find(parser.taxon());
                        if (it != taxon_map.end()) {
                                occurred.insert(parser.taxon());
                                sequences[it->second] = sequence_t<AC>(line, alphabet);
                        }
                        else {
                                std::cerr << boost::format("Warning: taxon `%s' not found in the phylogenetic tree.") % parser.taxon()
                                          << std::endl;
                        }
                }
                push_back(alignment_t<AC>(sequences, taxon_map, alphabet, verbose));
                if (verbose) {
                        std::cerr << boost::format("Finished parsing alignment %d...") % (base_t::size())
                                  << std::endl;
                }
        }
        // use the access operator from base class to read/write columns
        using base_t::operator[];
        alignment_t<AC> operator[](const range_t& range) const {
                index_t index(range.index()[1]);
                range_t tmp(index, range.length());
                return base_t::operator[](range.index()[0])[tmp];
        }

protected:
        // Methods for Initializations
        ////////////////////////////////////////////////////////////////////////
        taxon_map_t create_taxon_map(const pt_root_t& tree) {
                taxon_map_t taxon_map;

                for (pt_node_t::leaves_t::const_iterator it = tree.begin_leaves();
                     it != tree.end_leaves(); it++) {
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

// data structure optimized for computing likelihoods
////////////////////////////////////////////////////////////////////////////////
template <typename AC = alphabet_code_t>
class alignment_map_t : public std::map<std::vector<AC>, double>
{
        typedef std::map<std::vector<AC>, double> base_t;
public:
        explicit alignment_map_t(const alignment_t<AC>& alignment)
                : base_t() {
                for (typename alignment_t<AC>::const_iterator it = alignment.begin();
                     it != alignment.end(); it++) {
                        base_t::operator[](*it) += 1.0;
                }
        }
        explicit alignment_map_t(const alignment_set_t<AC>& alignment_set)
                : base_t() {
                for (typename alignment_set_t<AC>::const_iterator it = alignment_set.begin();
                     it != alignment_set.end(); it++) {
                        for (typename alignment_t<AC>::const_iterator is = it->begin();
                             is != it->end(); is++) {
                                base_t::operator[](*is) += 1.0;
                        }
                }
                size_t s = 0;
                for (typename alignment_set_t<AC>::const_iterator it = alignment_set.begin();
                     it != alignment_set.end(); it++) {
                        s += it->length();
                }
        }
};

// i/o streams
////////////////////////////////////////////////////////////////////////////////
template <typename AC>
class print_alignment_t {
public:
        print_alignment_t(const alignment_t<AC>& alignment)
                : _alignment(alignment)
                { }

        virtual const alignment_t<AC>& operator()() const {
                return _alignment;
        }
protected:
        const alignment_t<AC>& _alignment;
};

template <typename AC>
class print_alignment_pretty_t : public print_alignment_t<AC> {
public:
        print_alignment_pretty_t(const alignment_t<AC>& alignment)
                : print_alignment_t<AC>(alignment)
                { }
};
template <typename AC>
print_alignment_pretty_t<AC> print_alignment_pretty(alignment_t<AC> alignment)
{
        return print_alignment_pretty_t<AC>(alignment);
}

template <typename AC>
class print_alignment_fasta_t : public print_alignment_t<AC> {
public:
        print_alignment_fasta_t(const alignment_t<AC>& alignment)
                : print_alignment_t<AC>(alignment)
                { }
};
template <typename AC>
print_alignment_fasta_t<AC> print_alignment_fasta(alignment_t<AC> alignment)
{
        return print_alignment_fasta_t<AC>(alignment);
}

template <typename AC>
std::ostream& operator<< (std::ostream& o, const print_alignment_pretty_t<AC>& alignment_container)
{
        const alignment_t<AC>& alignment(alignment_container());

        for (typename alignment_t<AC>::taxon_map_t::const_iterator it = alignment.taxon_map().begin();
             it != alignment.taxon_map().end(); it++) {

                o << boost::format("%10s (%2d): ") % it->first % it->second
                  << alignment[it->first]
                  << std::endl;
        }
        return o;
}

template <typename AC>
std::ostream& operator<< (std::ostream& o, const print_alignment_fasta_t<AC>&  alignment_container)
{
        const alignment_t<AC>& alignment(alignment_container());

        for (typename alignment_t<AC>::taxon_map_t::const_iterator it = alignment.taxon_map().begin();
             it != alignment.taxon_map().end(); it++) {
                std::stringstream ss;
                ss << alignment[it->first];
                o << boost::format(">%s\n") % it->first
                  << split_string(ss.str(), 80)
                  << std::endl;
        }
        return o;
}

#endif /* __TFBAYES_ALIGNMENT_ALIGNMENT_HH__ */
