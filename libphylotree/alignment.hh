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
        alignment_t(const char* filename, pt_root_t* tree)
                : tree(tree) {

                FastaParser parser(filename);
                std::string sequence;

                tree->init(4);

                while ((sequence = parser.read_sequence()) != "") {
                        std::string taxon = token(parser.description()[0], '.')[0];
                        pt_node_t* node   = find_node(taxon, tree);
                        if (node) {
                                taxa_mapping.insert(std::pair<std::string, pt_node_t*>(taxon, node));
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
                        length = alignments[taxa_mapping.begin()->first].size();
                }
                else {
                        length = 0;
                }
        }
        pt_node_t* find_node(const std::string& name, pt_node_t* node) {
                if (node->name == name) {
                        return node;
                }
                else {
                        if (node->leaf()) {
                                return NULL;
                        }
                        pt_node_t* tmp = find_node(name, node->left);
                        if (tmp) {
                                return tmp;
                        }
                        return find_node(name, node->right);
                }
        }

        static
        const std::string strip(const std::string& str) {
                const std::string& whitespace = " \t\n";
                const size_t begin = str.find_first_not_of(whitespace);
                if (begin == std::string::npos) {
                        return "";
                }
                const size_t end   = str.find_last_not_of(whitespace);
                const size_t range = end - begin + 1;

                return str.substr(begin, range);
        }

        static
        std::vector<std::string> token(const std::string& str, char t) {
                std::string token;
                std::vector<std::string> tokens;
                std::istringstream iss(str);
                while (getline(iss, token, t)) {
                        tokens.push_back(strip(token));
                }
                return tokens;
        }

        // typedefs
        ////////////////////////////////////////////////////////////////////////
        typedef boost::unordered_map<std::string, nucleotide_sequence_t<CODE_TYPE> > alignment_type;
        typedef std::set<std::pair<std::string, pt_node_t*> > taxa_mapping_type;

        // Iterator
        ////////////////////////////////////////////////////////////////////////
        class iterator
        {
        public:
                iterator(pt_root_t* tree, size_t i, size_t length,
                         const alignment_type& alignments,
                         const taxa_mapping_type& taxa_mapping)
                        : tree(tree), i(i), length(length), alignments(alignments),
                          taxa_mapping(taxa_mapping) { 
                        //
                        insert_observables();
                }

                void insert_observables() {
                        if (i >= length) {
                                return;
                        }
                        for (taxa_mapping_type::const_iterator it = taxa_mapping.begin();
                             it != taxa_mapping.end(); it++) {
                                it->second->x = alignments.find(it->first)->second[i];
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
                        insert_observables();
                        return *this;
                }
                pt_root_t* operator->() {
                        return tree;
                }
                pt_root_t& operator*() {
                        return *operator->();
                }
        private:
                pt_root_t* tree;
                size_t i;
                size_t length;
                const alignment_type& alignments;
                const taxa_mapping_type& taxa_mapping;
        };
        iterator begin() {
                return iterator(tree, 0, length, alignments, taxa_mapping);
        }
        iterator end() {
                return iterator(tree, length, length, alignments, taxa_mapping);
        }

        pt_root_t* tree;
        size_t length;
private:
        alignment_type alignments;
        taxa_mapping_type taxa_mapping;
};

#endif /* ALIGNMENT_HH */
