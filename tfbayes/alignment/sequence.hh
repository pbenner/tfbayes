/* Copyright (C) 2011-2013 Philipp Benner
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

#ifndef __TFBAYES_ALIGNMENT_SEQUENCE_HH__
#define __TFBAYES_ALIGNMENT_SEQUENCE_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <functional>
#include <string>
#include <vector>

#include <boost/bind.hpp>            /* bind */
#include <boost/range/adaptors.hpp>  /* transformed */
#include <boost/range/algorithm.hpp> /* copy */
#include <boost/foreach.hpp>

#include <tfbayes/uipac/alphabet.hh>

template <typename CODE_TYPE = alphabet_code_t>
class sequence_t : public std::vector<CODE_TYPE>
{
public:
        sequence_t(const alphabet_t& alphabet = nucleotide_alphabet_t())
                : std::vector<CODE_TYPE>(),
                  alphabet(alphabet)
                { }
        sequence_t(const size_t length,
                   const alphabet_t& alphabet = nucleotide_alphabet_t())
                : std::vector<CODE_TYPE>(length, alphabet.code('N')),
                  alphabet(alphabet)
                { }
        sequence_t(const size_t length, CODE_TYPE init,
                   const alphabet_t& alphabet = nucleotide_alphabet_t())
                : std::vector<CODE_TYPE>(length, alphabet.code(init)),
                  alphabet(alphabet)
                { }
        sequence_t(const std::string& sequence,
                   const alphabet_t& alphabet = nucleotide_alphabet_t())
                : std::vector<CODE_TYPE>(sequence.size()),
                  alphabet(alphabet) {
                boost::copy(
                        sequence | boost::adaptors::transformed(boost::bind(&alphabet_t::code, &alphabet, _1)), begin());
        }
        sequence_t(const sequence_t<CODE_TYPE>& ns,
                   const alphabet_t& alphabet = nucleotide_alphabet_t())
                : std::vector<CODE_TYPE>(ns),
                  alphabet(alphabet) {
        }
        sequence_t(const std::vector<CODE_TYPE>& ns,
                   const alphabet_t& alphabet = nucleotide_alphabet_t())
                : std::vector<CODE_TYPE>(ns),
                  alphabet(alphabet) {
        }
        sequence_t(const sequence_t<CODE_TYPE>& ns, const std::string& sequence,
                   const alphabet_t& alphabet = nucleotide_alphabet_t())
                : std::vector<CODE_TYPE>(ns),
                  alphabet(alphabet) {
                for (size_t i = 0; i < sequence.length(); i++) {
                        this->push_back(alphabet.code(sequence[i]));
                }
        }
        using std::vector<CODE_TYPE>::begin;
        using std::vector<CODE_TYPE>::end;

        friend std::ostream& operator<< (std::ostream& o, const sequence_t<CODE_TYPE>& sequence) {
                for (typename sequence_t<CODE_TYPE>::const_iterator it = sequence.begin();
                     it != sequence.end(); it++) {
                        o << (*it == sequence.alphabet.code('N') ? 'N' : sequence.alphabet.decode(*it));
                }
                return o;
        }
protected:
        alphabet_t alphabet;
};

#endif /* __TFBAYES_ALIGNMENT_SEQUENCE_HH__ */
