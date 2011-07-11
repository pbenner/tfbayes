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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <string.h>

#include <code.hh>

using namespace std;

bool is_nucleotide(char S) {
        if (S == 'A' || S == 'C' || S == 'G' || S == 'T' ||
            S == 'a' || S == 'c' || S == 'g' || S == 't') {
                return true;
        }
        else {
                return false;
        }
}

char code_nucleotide(char a)
{
        switch (a) {
        case 'A':
        case 'a':
                return 0;
        break;
        case 'C':
        case 'c':
                return 1;
        break;
        case 'G':
        case 'g':
                return 2;
        break;
        case 'T':
        case 't':
                return 3;
        break;
        default:
                return -1;
        }
}

char decode_nucleotide(char a)
{
        switch (a) {
        case 0:
                return 'a';
        break;
        case 1:
                return 'c';
        break;
        case 2:
                return 'g';
        break;
        case 3:
                return 't';
        break;
        }
        std::cerr << "decode_nucleotide(): internal error.\n"
                  << std::endl;
        exit(EXIT_FAILURE);
}

vector<char> code_nucleotide_sequence(const char* sequence) throw(InvalidNucleotide)
{
        size_t n = strlen(sequence);
        vector<char> coded_sequence(n, 0);

        for (size_t i = 0; i < n; i++) {
                char c = code_nucleotide(sequence[i]);
                if (c == -1) {
                        throw InvalidNucleotide("code_nucleotide_sequence()", sequence[i], sequence);
                }
                coded_sequence[i] = c;
        }
        return coded_sequence;
}
