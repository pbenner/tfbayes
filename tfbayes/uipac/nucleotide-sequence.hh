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

#ifndef NUCLEOTIDE_SEQUENCE_HH
#define NUCLEOTIDE_SEQUENCE_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <string>

#include <tfbayes/uipac/code.hh>

template <typename CODE_TYPE>
class nucleotide_sequence_t : public std::vector<CODE_TYPE>
{
public:
        nucleotide_sequence_t()
                : std::vector<CODE_TYPE>() { }
        nucleotide_sequence_t(const size_t length)
                : std::vector<CODE_TYPE>(length, code_nucleotide<CODE_TYPE>('-')) { }
        nucleotide_sequence_t(const size_t length, CODE_TYPE init)
                : std::vector<CODE_TYPE>(length, init) { }
        nucleotide_sequence_t(const std::string& sequence)
                : std::vector<CODE_TYPE>() {
                for (size_t i = 0; i < sequence.length(); i++) {
                        this->push_back(code_nucleotide<CODE_TYPE>(sequence[i]));
                }
        }
};

#endif /* NUCLEOTIDE_SEQUENCE_HH */
