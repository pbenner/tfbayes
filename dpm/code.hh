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

#ifndef CODE_HH
#define CODE_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <string>
#include <vector>

class InvalidNucleotide: public std::exception
{
public:
        InvalidNucleotide(const char* fname, char nucleotide, const std::string sequence)
                : fname(fname), _nucleotide(nucleotide), _sequence(sequence)
                { }
        ~InvalidNucleotide() throw() {};

        virtual const char* what() const throw() {
                return (std::string(fname) + std::string(": found invalid nucleotide")).c_str();
        }
        const std::string& sequence() const {
                return _sequence;
        }
        char nucleotide() const {
                return _nucleotide;
        }

private:
        const char* fname;
        const char _nucleotide;
        const std::string _sequence;
};

bool is_nucleotide(char a);
char code_nucleotide(char a);
char decode_nucleotide(char a);
std::vector<char> code_nucleotide_sequence(const std::string& sequence) throw(InvalidNucleotide);

#endif /* CODE_HH */
