/* Copyright (C) 2015 Philipp Benner
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

#include <iostream>
#include <cerrno>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <boost/format.hpp>

#include <tfbayes/config/partition-parser.hh>

using namespace std;

int yylex_init   (yyscan_t* scanner);
int yylex_destroy(yyscan_t  scanner);
void yylex_set_input(yyscan_t scanner, FILE* file);

dpm_partition_list_t parse_partition_list(FILE * file)
{
        context_t context;

        // initialize lexer
        yylex_init(&context.scanner);

        // set input file if necessary
        if (file != NULL) {
                yylex_set_input(context.scanner, file);
        }

        // parse input
        if (yyparse(&context)) {
                // error parsing
                return context.partition_list;
        }
 
        // free lexer memory
        yylex_destroy(context.scanner);

        return context.partition_list;
}

dpm_partition_list_t parse_partition_list(const string& filename)
{
        FILE* yyin = NULL;
        if (filename != "") {
                yyin = fopen(filename.c_str(), "r");
                if (yyin == NULL) {
                        cerr << boost::format("Could not open file `%s': %s") % filename % strerror(errno)
                             << endl;
                        exit(EXIT_FAILURE);
                }
        }
        dpm_partition_list_t partition_list= parse_partition_list(yyin);
        if (yyin) fclose(yyin);

        return partition_list;
}
