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

        if (locp) {
                std::cerr << "locp: " << locp << std::endl;
                std::cerr << "text: " << yyget_text(context->scanner) << std::endl;

                int i = locp->first_line;
                int j = locp->first_column;
                if (i)  {
                        std::cerr << "ERROR 1" << std::endl;
                }
                if (j)  {
                        std::cerr << "ERROR 2" << std::endl;
                }
        }

        ss << boost::format("parsing error at line %d column %d near `%s': %s\n")
		% locp->first_line % locp->first_column
		% yyget_text(context->scanner) % std::string(err);

        std::cerr << ss.str() << std::endl;

        throw std::runtime_error(ss.str());
}
