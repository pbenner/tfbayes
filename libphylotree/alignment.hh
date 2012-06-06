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
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/unordered_map.hpp>
#include <set>
#include <string>

#include <phylotree.hh>
#include <tfbayes/fasta.hh>
#include <tfbayes/strtools.hh>

template <typename CODE_TYPE>
CODE_TYPE code_nucleotide(char a)
{
        switch (a) {
        case 'A':
        case 'a':
                return 1;
        break;
        case 'C':
        case 'c':
                return 3;
        break;
        case 'G':
        case 'g':
                return 0;
        break;
        case 'T':
        case 't':
                return 2;
        break;
        case 'N':
        case 'n':
        case '-':
                return 4;
        default:
                break;
        }
        return -1;
}

template <typename CODE_TYPE>
class nucleotide_sequence_t : public std::vector<CODE_TYPE>
{
public:
        nucleotide_sequence_t()
                : std::vector<CODE_TYPE>() { }
        nucleotide_sequence_t(const std::string& sequence)
                : std::vector<CODE_TYPE>() {
                for (size_t i = 0; i < sequence.length(); i++) {
                        push_back(code_nucleotide<CODE_TYPE>(sequence[i]));
                }
        }
};

template <typename CODE_TYPE>
class alignment_t {
public:
        alignment_t(const char* filename, pt_root_t* tree) {

                FastaParser parser(filename);
                std::string sequence;

                while ((sequence = parser.read_sequence()) != "") {
                        std::string taxon = token(parser.description()[0], '.')[0];
                        pt_node_t* node   = tree->map[taxon];
                        if (node) {
                                taxa.insert(taxon);
                                alignments[taxon] = nucleotide_sequence_t<CODE_TYPE>(sequence);
                        }
                        else {
                                std::cerr << "Warning: taxon `"
                                          << taxon
                                          << "' not found in the phylogenetic tree."
                                          << std::endl;
                        }
                }
                if (alignments.size() > 0) {
                        length = alignments[*taxa.begin()].size();
                }
                else {
                        length = 0;
                }
        }

        // typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_map<std::string, nucleotide_sequence_t<CODE_TYPE> > alignment_type;
        typedef std::set<std::string> taxa_type;

        // Iterator
        ////////////////////////////////////////////////////////////////////////
        class iterator
        {
        public:
                iterator(size_t i, size_t length,
                         const alignment_type& alignments,
                         const taxa_type& taxa)
                        : i(i), length(length), alignments(alignments),
                          taxa(taxa) { }

                void apply(pt_root_t* tree) {
                        if (i >= length) {
                                return;
                        }
                        for (taxa_type::const_iterator it = taxa.begin();
                             it != taxa.end(); it++) {
                                tree->map[*it]->x = alignments.find(*it)->second[i];
                        }
                }

                bool operator==(const iterator& it) {
                        return i == it.i;
                }
                bool operator!=(const iterator& it) {
                        return i != it.i;
                }
                iterator& operator++(int _) {
                        i++;
                        return *this;
                }
                // pt_root_t* operator->() {
                //         return tree;
                // }
                // pt_root_t& operator*() {
                //         return *operator->();
                // }
        private:
                size_t i;
                size_t length;
                const alignment_type& alignments;
                const taxa_type& taxa;
        };
        iterator begin() {
                return iterator(0, length, alignments, taxa);
        }
        iterator end() {
                return iterator(length, length, alignments, taxa);
        }

        // Reverse Iterator
        ////////////////////////////////////////////////////////////////////////
        class reverse_iterator
        {
        public:
                reverse_iterator(size_t i, size_t length,
                         const alignment_type& alignments,
                         const taxa_type& taxa)
                        : i(i), length(length), alignments(alignments),
                          taxa(taxa) { }

                void apply(pt_root_t* tree) {
                        if (i >= length) {
                                return;
                        }
                        for (taxa_type::const_iterator it = taxa.begin();
                             it != taxa.end(); it++) {
                                tree->map[*it]->x = alignments.find(*it)->second[i];
                        }
                }

                bool operator==(const reverse_iterator& it) {
                        return i == it.i;
                }
                bool operator!=(const reverse_iterator& it) {
                        return i != it.i;
                }
                reverse_iterator& operator++(int _) {
                        i++;
                        return *this;
                }
                // pt_root_t* operator->() {
                //         return tree;
                // }
                // pt_root_t& operator*() {
                //         return *operator->();
                // }
        private:
                size_t i;
                size_t length;
                const alignment_type& alignments;
                const taxa_type& taxa;
        };
        reverse_iterator rbegin() {
                return reverse_iterator(0, length, alignments, taxa);
        }
        reverse_iterator rend() {
                return reverse_iterator(length, length, alignments, taxa);
        }

        size_t length;
private:
        alignment_type alignments;
        taxa_type taxa;
};

#endif /* ALIGNMENT_HH */
