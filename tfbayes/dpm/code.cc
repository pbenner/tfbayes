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

bool is_nucleotide_or_masked(char S) {
        if (S == 'A' || S == 'C' || S == 'G' || S == 'T' || S == 'N' ||
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
                return 4;
        default:
                break;
        }
        return -1;
}

char decode_nucleotide(char a)
{
        switch (a) {
        case 1:
                return 'A';
        break;
        case 3:
                return 'C';
        break;
        case 0:
                return 'G';
        break;
        case 2:
                return 'T';
        break;
        case 4:
                return 'N';
        default:
                break;
        }
        std::cerr << "decode_nucleotide(): internal error.\n"
                  << std::endl;
        exit(EXIT_FAILURE);
}

char complement_nucleotide(char a)
{
        switch (a) {
        case 'A':
                return 'T';
        case 'a':
                return 't';
        case 'C':
                return 'G';
        case 'c':
                return 'g';
        case 'G':
                return 'C';
        case 'g':
                return 'c';
        case 'T':
                return 'A';
        case 't':
                return 'a';
        case 'N':
                return 'N';
        default:
                break;
        }
        return -1;
}

vector<short> code_nucleotide_sequence(const string& sequence) throw(InvalidNucleotide)
{
        size_t n = sequence.size();
        vector<short> coded_sequence(n, 0);

        for (size_t i = 0; i < n; i++) {
                char c = code_nucleotide(sequence[i]);
                if (c == -1) {
                        throw InvalidNucleotide("code_nucleotide_sequence()", sequence[i], sequence);
                }
                coded_sequence[i] = c;
        }
        return coded_sequence;
}

vector<vector<short> > code_sequences(const vector<string>& sequences) {
        vector<vector<short> > sequences_coded;

        for(size_t i = 0; i < sequences.size(); i++) {
                try {
                        sequences_coded.push_back(code_nucleotide_sequence(sequences[i]));
                }
                catch (const InvalidNucleotide& i) {
                        cout << "Exception:    " << i.what()     << " \'" << i.nucleotide() << "\'"
                             << " in sequence: " << i.sequence() << endl;
                        exit(EXIT_FAILURE);
                }
        }
        return sequences_coded;
}

vector<string> complement(const vector<string>& sequences) {
        vector<string> sequences_reversed;

        for(size_t i = 0; i < sequences.size(); i++) {
                sequences_reversed.push_back(string(sequences[i].size(), '\0'));
                for(size_t j = 0; j < sequences[i].size(); j++) {
                        sequences_reversed[i][j] = complement_nucleotide(sequences[i][j]);
                }
        }
        return sequences_reversed;
}
