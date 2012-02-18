%{
#include <stdio.h>
#include "phylotree-parser.h"
size_t line_count;
extern int yylval;
%}

WHITESPACE [\ \t]
NEWLINE [\n]

%%
{WHITESPACE}+              { }
{NEWLINE}                  { line_count++; }
root                       { return ROOT; }
node                       { return NODE; }
leaf                       { return LEAF; }
[a-zA-Z][a-zA-Z0-9]*       { return ID; }
\)                         { return RPAREN; }
\(                         { return LPAREN; }
[0-9]+                     { yylval = atoi(yytext); return NUM; }
\"[^\"\n]*\"               { return STR; }
.                          { return yytext[0]; }
%%
int yywrap(void) {
        return 1;
}
