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

FastaParser::FastaParser(const string& filename)
        : file(filename.c_str()),
          stream(file)
{
        if (!file.good())
        {
                cerr << "Could not open `"
                     << filename
                     << "'."
                     << endl;
        }
        init_parser();
}

FastaParser::FastaParser(const istream& stream)
        : stream(file)
{
        init_parser();
}

void
FastaParser::init_parser()
{
        string tmp;

        while (stream) {
                getline(stream, tmp);
                string line = strip(tmp);
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
FastaParser::description() const
{
        return _description;
}

string
FastaParser::taxon() const
{
        if (description().size() == 0)
                return "";
        else {
                return token(description()[0], '.')[0];
        }
}

string
FastaParser::operator()()
{
        string tmp;
        string sequence;

        if (prev_line == "" || !stream) {
                return "";
        }
        if (prev_line[0] != '>') {
                cerr << "Invalid fasta file." << endl;
                exit(EXIT_FAILURE);
        }
        prev_line.erase(0, 1);
        _description = token(prev_line, '|');

        while (true) {
                getline(stream, tmp);
                string line = strip(tmp);
                if (line == "" && stream) {
                        continue;
                }
                if (tmp[0] == '>' || !stream) {
                        prev_line = line;
                        break;
                }
                sequence.append(line);
        }
        return sequence;
}

FastaParser::operator bool() const
{
        return stream.good();
}
