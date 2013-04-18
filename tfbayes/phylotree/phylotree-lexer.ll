
/* make the scanner thread-safe (no global variables) */
%option reentrant
/* tell flex that the scanner is combined with a bison parser */
%option bison-bridge
%option bison-locations
%option yylineno

%{
#include <stdio.h>
#include <tfbayes/phylotree/phylotree-parser.h>

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
{NEWLINE}                  { }
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
