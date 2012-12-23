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

#include <vector>
#include <string>

#include <math.h>

#include <tfbayes/utility/abysmal-stack.hh>
#include <tfbayes/dpm/code.hh>

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

class context_t : public std::vector<int>
{
public:
        context_t(size_t alphabet_size, const AbysmalStack<short>& stack) {
                for (size_t context = 0; context < stack.depth(); context++) {
                        if (!stack.clogged(context)) {
                                // compute counts position
                                const size_t offset   = counts_offset(alphabet_size, context);
                                const size_t position = counts_position(alphabet_size, context, stack);
                                push_back(offset+position);
                        }
                        else {
                                // no context available
                                push_back(-1);
                        }
                }
        }

        static size_t counts_size(size_t alphabet_size, size_t max_context) {
                return counts_offset(alphabet_size, max_context+1);
        }

        static size_t counts_offset(size_t alphabet_size, size_t context) {
                size_t offset = 0;

                for (size_t i = 0; i < context+1; i++) {
                        offset += pow(alphabet_size, i);
                }

                return offset - 1;
        }

protected:
        static size_t counts_position(size_t alphabet_size, size_t context, const AbysmalStack<short>& stack) {
                size_t position = 0;

                for (size_t i = 0; i < context; i++) {
                        position += stack[i+stack.depth()-1-context]*pow(alphabet_size, i+1);
                }

                return position + stack.top();
        }
};


class seq_context_t : public std::vector<context_t>
{
public:
        seq_context_t(const nucleotide_sequence_t& sequence,
                      size_t context, size_t alphabet_size) {
                AbysmalStack<short> stack(context+1);

                for (size_t i = 0; i < sequence.size(); i++) {
                        if (sequence[i] == (short)alphabet_size) {
                                stack.push_invalid(sequence[i]);
                        }
                        else {
                                stack.push(sequence[i]);
                        }
                        push_back(context_t(alphabet_size, stack));
                }
        }
};

#endif /* NUCLEOTIDE_SEQUENCE_HH */
