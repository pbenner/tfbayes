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

#ifndef FASTA_HH
#define FASTA_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

class FastaParser {
public:
        FastaParser(const std::string& filename);
        FastaParser(const std::istream& stream);

        /* read the next sequence */
        std::string operator()();
        /* return true if not at end of file */
        operator bool() const;
        /* return the description of the sequence */
        const std::vector<std::string>& description() const;

protected:
        void init_parser();

        std::ifstream file;
        std::istream& stream;
        std::string prev_line;
        std::vector<std::string> _description;
};

#endif /* FASTA_HH */
