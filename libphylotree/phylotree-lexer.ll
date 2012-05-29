%{
#include <stdio.h>
#include "phylotree-parser.h"
size_t line_count;
extern int yylval;
%}

COLON :
COMMA ,
SEMICOLON ;
DIGIT [0-9]
NEWLINE [\n]
WHITESPACE [\ \t]

%%
{WHITESPACE}+              { }
{NEWLINE}                  { line_count++; }
{COLON}                    { return COLON; }
{COMMA}                    { return COMMA; }
{SEMICOLON}                { return SEMICOLON; }
[a-zA-Z][a-zA-Z0-9]*       { return NAME; }
\)                         { return RPAREN; }
\(                         { return LPAREN; }
-?{DIGIT}+("."{DIGIT}*)?   { return FLOAT; }
.                          { return yytext[0]; }
%%
int yywrap(void) {
        return 1;
}
