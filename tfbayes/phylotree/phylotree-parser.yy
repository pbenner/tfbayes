
// lexer definitions
////////////////////////////////////////////////////////////////////////////////

%code requires {
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/phylotree/phylotree-parsetree.hh>
#define YYSTYPE pt_parsetree_t *
#define YYLEX_PARAM context->scanner
}

%code {

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tfbayes/phylotree/phylotree-parser.h>

typedef void* yyscan_t;
char *yyget_text (yyscan_t scanner);

int yylex(YYSTYPE * yylval_param, YYLTYPE * yylloc_param, yyscan_t scanner);

int yyerror(YYLTYPE* locp, context_t* context, const char* err) {
        fprintf(stderr, "parsing error at line %d colum %d near `%s': %s\n",
		locp->first_line+1, locp->first_column,
		yyget_text(context->scanner), err);
        exit(EXIT_FAILURE);
}

}

// bison options
////////////////////////////////////////////////////////////////////////////////

%pure-parser
%locations
%defines
%error-verbose
%parse-param {context_t* context}
%lex-param {void * scanner}

// token definitions
////////////////////////////////////////////////////////////////////////////////

%token COLON COMMA SEMICOLON LPAREN RPAREN NAME FLOAT

// grammar
////////////////////////////////////////////////////////////////////////////////
%%
start:
      tree_list root SEMICOLON
      { context->pt_parsetree = new pt_parsetree_t(TREE_N, 2, NULL, $1, $2); }
    | root SEMICOLON
      { context->pt_parsetree = new pt_parsetree_t(TREE_N, 1, NULL, $1); }
    ;
tree_list:
      tree_list root SEMICOLON
      { $$ = new pt_parsetree_t(TREE_N, 2, NULL, $1, $2); }
    | root SEMICOLON
      { $$ = new pt_parsetree_t(TREE_N, 1, NULL, $1); }
    ;
root:  LPAREN node COMMA node COMMA node RPAREN
      { $$ = new pt_parsetree_t(ROOT_N, 3, NULL, $2, $4, $6); }
    | LPAREN node COMMA node RPAREN
      { $$ = new pt_parsetree_t(ROOT_N, 2, NULL, $2, $4); }
    ;
node: LPAREN node COMMA node RPAREN COLON distance
      { $$ = new pt_parsetree_t(NODE_N, 3, NULL, $2, $4, $7); }
    | name COLON distance
      { $$ = new pt_parsetree_t(LEAF_N, 2, NULL, $1, $3); }
    ;
name: NAME { $$ = new pt_parsetree_t(NAME_N, 0, strdup(yyget_text(context->scanner))); }
    ;
distance:
      FLOAT
      {
        $$ = new pt_parsetree_t(DISTANCE_N, 0, calloc(1, sizeof(double)));
        *((double *)$$->data) = atof(yyget_text(context->scanner));
      }
    ;
%%
