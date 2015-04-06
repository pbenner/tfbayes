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
#include <sstream>
#include <stdexcept>
#include <cerrno>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <boost/format.hpp>

using namespace std;

// tools
////////////////////////////////////////////////////////////////////////////////

int yywrap(yyscan_t scanner) {
        return 1;
}

void yylex_set_input(yyscan_t scanner, FILE* file) {
	struct yyguts_t * yyg = (struct yyguts_t*)scanner;
	yyg->yyin_r = file;
}

// interface
////////////////////////////////////////////////////////////////////////////////

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
        FILE* fp = NULL;
        if (filename != "") {
                fp = fopen(filename.c_str(), "r");
                if (fp == NULL) {
                        cerr << boost::format("Could not open file `%s': %s") % filename % strerror(errno)
                             << endl;
                        exit(EXIT_FAILURE);
                }
        }
        dpm_partition_list_t partition_list = parse_partition_list(fp);
        if (fp) fclose(fp);

        return partition_list;
}

dpm_partition_list_t parse_partition_list_str(const string& str)
{
        context_t context;

        // initialize lexer
        yylex_init(&context.scanner);

        // read buffer
        YY_BUFFER_STATE buf = yy_scan_string(str.c_str(), context.scanner);

        // it seems that the line and column numbers have to be
        // initialized manually
        buf->yy_bs_lineno = 1;
        buf->yy_bs_column = 0;

        // parse input
        yyparse(&context);

        // delete string buffer
        yy_delete_buffer(buf, context.scanner);
 
        // free lexer memory
        yylex_destroy(context.scanner);

        return context.partition_list;
}

dpm_partition_t parse_partition_str(const string& str)
{
        dpm_partition_list_t partition_list = parse_partition_list_str(str);

        if (partition_list.size() != 1) {
                stringstream ss;
                ss << boost::format("Error: Found %d partitions but expected one.")
                        % partition_list.size();
                throw runtime_error(ss.str());
        }
        return *partition_list.begin();
}
