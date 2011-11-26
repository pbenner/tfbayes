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
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <string>

#include <abysmal-stack.hh>
#include <code.hh>

#include <parsmm/static_pars_tree.h>

class nucleotide_sequence_t : public std::vector<short>
{
public:
        nucleotide_sequence_t(const std::string& sequence)
                : _sequence(sequence) {
                for (size_t i = 0; i < sequence.length(); i++) {
                        push_back(code_nucleotide(sequence[i]));
                }
        }

        const std::string& toString() {
                return _sequence;
        }

private:
        const std::string _sequence;
};

class nucleotide_context_t : public std::vector<short>
{
public:
        nucleotide_context_t(const nucleotide_sequence_t& sequence,
                             size_t depth, size_t alphabet_size) {
                AbysmalStack<count_t> stack(depth+1);
                size_t position;

                for (size_t i = 0; i < sequence.size(); i++) {
                        stack.push(sequence[i]);
                        if (i >= depth) {
                                position = pow(alphabet_size, depth+1) - 1;
                                for (size_t i = 0; i < depth+1; i++) {
                                        position -= stack[i]*pow(alphabet_size, i);
                                }
                                std::cout << stack << " -> " << position << std::endl;
                                push_back(position);
                        }
                        else {
                                push_back(-1);
                        }
                }
        }
};

#endif /* NUCLEOTIDE_SEQUENCE_HH */
