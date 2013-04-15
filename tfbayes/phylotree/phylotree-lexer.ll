
/* make the scanner thread-safe (no global variables) */
%option reentrant
/* tell flex that the scanner is combined with a bison parser */
%option bison-bridge
%option bison-locations

%{
#include <stdio.h>
#include <tfbayes/phylotree/phylotree-parser.h>
size_t line_count;
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
int yywrap(yyscan_t scanner) {
        return 1;
}
void yylex_set_input(yyscan_t scanner, FILE* file) {
	struct yyguts_t * yyg = (struct yyguts_t*)scanner;
	yyg->yyin_r = file;
}
