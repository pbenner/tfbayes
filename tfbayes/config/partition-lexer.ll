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

/* make the scanner thread-safe (no global variables) */
%option reentrant
/* tell flex that the scanner is combined with a bison parser */
%option bison-bridge
%option bison-locations
%option yylineno
%option nounput

%{
#include <tfbayes/config/partition-parser.hh>

#define YY_USER_ACTION \
	yylloc->first_line   = yylloc->last_line = yylineno; \
	yylloc->first_column = yycolumn; \
	yylloc->last_column  = yycolumn + yyleng - 1; \
	yycolumn            += yyleng;

%}

COLON :
COMMA ,
SEMICOLON ;
DIGIT [0-9]
NEWLINE [\n]
WHITESPACE [\ \t]

%%
{WHITESPACE}+              { }
{NEWLINE}                  { return NEWLINE;}
{COLON}                    { return COLON; }
{COMMA}                    { return COMMA; }
{SEMICOLON}                { return SEMICOLON; }
[a-zA-Z_][a-zA-Z0-9_-]*    { return STRING; }
\)                         { return RPAREN; }
\(                         { return LPAREN; }
\}                         { return RCPAREN; }
\{                         { return LCPAREN; }
{DIGIT}+                   { return INTEGER; }
!                          { return EXCLAMATION; }
-                          { return DASH; }
.                          { return yytext[0]; }
%%
#include <tfbayes/config/partition-lexer-code.hh>
