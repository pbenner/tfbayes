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

#include <new>
#include <inttypes.h>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include <boost/format.hpp>

#include <tfbayes/config/partition-parser.hh>

char *yyget_text (yyscan_t scanner);

int yylex(YYSTYPE * yylval_param, YYLTYPE * yylloc_param, yyscan_t scanner);
int yylex(YYSTYPE * yylval_param, YYLTYPE * yylloc_param, context_t* context) {
    return yylex(yylval_param, yylloc_param, context->scanner);
}

int yyerror(YYLTYPE* locp, context_t* context, const char* err) {
        std::stringstream ss;

        if (*yyget_text(context->scanner)) {
                ss << boost::format("parsing error at line %d column %d near `%s': %s\n")
                        % locp->first_line % locp->first_column
                        % yyget_text(context->scanner) % std::string(err);
        }
        else {
                ss << boost::format("parsing error at line %d column %d: %s\n")
                        % locp->first_line % locp->first_column
                        % std::string(err);
        }

        throw std::runtime_error(ss.str());
}
