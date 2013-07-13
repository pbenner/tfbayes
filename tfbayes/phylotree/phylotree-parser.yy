/* Copyright (C) 2013 Philipp Benner
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
		locp->first_line, locp->first_column,
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
      tree_list tree SEMICOLON
      { context->pt_parsetree = new pt_parsetree_t(TREE_LIST_N, 2, NULL, $1, $2); }
    | tree SEMICOLON
      { context->pt_parsetree = new pt_parsetree_t(TREE_LIST_N, 1, NULL, $1); }
    ;
tree_list:
      tree_list tree SEMICOLON
      { $$ = new pt_parsetree_t(TREE_LIST_N, 2, NULL, $1, $2); }
    | tree SEMICOLON
      { $$ = new pt_parsetree_t(TREE_LIST_N, 1, NULL, $1); }
    ;
tree: LPAREN node_list COMMA outgroup RPAREN
      { $$ = new pt_parsetree_t(TREE_N, 2, NULL, $2, $4); }
    ;
outgroup:
      name COLON distance
      { $$ = new pt_parsetree_t(LEAF_N, 2, NULL, $1, $3); }
node_list:
      node_list COMMA node
      { $$ = new pt_parsetree_t(NODE_LIST_N, 2, NULL, $1, $3); }
    | node
      { $$ = new pt_parsetree_t(NODE_LIST_N, 1, NULL, $1); }
    ;
node: name COLON distance
      { $$ = new pt_parsetree_t(LEAF_N, 2, NULL, $1, $3); }
    | LPAREN node_list RPAREN COLON distance
      { $$ = new pt_parsetree_t(NODE_N, 2, NULL, $2, $5); }
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
