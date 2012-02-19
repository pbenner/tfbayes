%{
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

#include <phylotree-parsetree.hh>

#define YYSTYPE pt_parsetree_t *

int yylex (void);
void yyerror (const char *);

pt_parsetree_t* root;
extern char *yytext;

%}

%token COLON COMMA SEMICOLON LPAREN RPAREN NAME FLOAT

%%
root: LPAREN node COMMA node RPAREN SEMICOLON
      { root = new pt_parsetree_t(ROOT_N, 2, NULL, $2, $4);; }
    | LPAREN node COMMA node RPAREN
      { root = new pt_parsetree_t(ROOT_N, 2, NULL, $2, $4);; }
    ;
node: LPAREN node COMMA node RPAREN COLON distance
      { $$ = new pt_parsetree_t(NODE_N, 3, NULL, $2, $4, $7); }
    | name COLON distance
      { $$ = new pt_parsetree_t(LEAF_N, 2, NULL, $1, $3); }
    ;
name: NAME { $$ = new pt_parsetree_t(NAME_N, 0, strdup(yytext)); }
    ;
distance:
      FLOAT
      {
        $$ = new pt_parsetree_t(DISTANCE_N, 0, calloc(1, sizeof(double)));
        *((double *)$$->data) = atof(yytext);
      }
    ;
%%
