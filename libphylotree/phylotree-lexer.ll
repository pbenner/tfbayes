%{
#include <stdio.h>
#include "phylotree-parser.h"
size_t line_count;
extern int yylval;
%}
%%
\n                         { line_count++; }
[0-9]+                     { yylval = atoi(yytext); return NUM; }
\"[^\"\n]*\"               { return STR; }
[a-zA-Z][a-zA-Z0-9]*       { return ID; }
\)                         { return RPAREN; }
\(                         { return LPAREN; }
.
%%
int yywrap(void) {
        return 1;
}
