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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>

#include <tfbayes/fasta/fasta.hh>
#include <tfbayes/utility/strtools.hh>

using namespace std;

FastaParser::FastaParser(const string& file_name)
        : file(file_name.c_str())
{
        string _line;

        if (!file.good())
        {
                cerr << "Could not open `"
                     << file_name
                     << "'."
                     << endl;
        }

        while (file) {
                getline(file, _line);
                string line = strip(_line);
                if (line == "" || line[0] == '>') {
                        prev_line = line;
                        break;
                }
                else {
                        cerr << "Invalid fasta file." << endl;
                        exit(EXIT_FAILURE);
                }
        }
}

const std::vector<std::string>&
FastaParser::description()
{
        return _description;
}

string
FastaParser::read_sequence()
{
        string _line;
        string sequence;

        if (prev_line == "" || !file) {
                return "";
        }
        if (prev_line[0] != '>') {
                cerr << "Invalid fasta file." << endl;
                exit(EXIT_FAILURE);
        }
        prev_line.erase(0, 1);
        _description = token(prev_line, '|');

        while (true) {
                getline(file, _line);
                string line = strip(_line);
                if (line == "" || _line[0] == '>' || !file) {
                        prev_line = line;
                        break;
                }
                sequence.append(line);
        }
        return sequence;
}
