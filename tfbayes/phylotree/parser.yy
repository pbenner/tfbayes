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

// general options
////////////////////////////////////////////////////////////////////////////////

%pure-parser
%locations
%defines
%error-verbose
%lex-param   {context_t* context}
%parse-param {context_t* context}

// lexer definitions
////////////////////////////////////////////////////////////////////////////////

%code requires {
#include <tfbayes/phylotree/parser-code-requires.hh>
}

%code {
#include <tfbayes/phylotree/parser-code.hh>
}

// token definitions
////////////////////////////////////////////////////////////////////////////////

%token COLON COMMA SEMICOLON LPAREN RPAREN NAME FLOAT

// grammar
////////////////////////////////////////////////////////////////////////////////
%%
start:
      tree_list
      { context->pt_parsetree = $1; }
    ;
tree_list:
      tree_list tree SEMICOLON
      { allocate($$, pt_parsetree_t(TREE_LIST_N, 2, NULL, $1, $2)); }
    | tree SEMICOLON
      { allocate($$, pt_parsetree_t(TREE_LIST_N, 1, NULL, $1)); }
    ;
tree: LPAREN node_list RPAREN
      { allocate($$, pt_parsetree_t(TREE_N, 1, NULL, $2)); }
    ;
node_list:
      node_list COMMA node
      { allocate($$, pt_parsetree_t(NODE_LIST_N, 2, NULL, $1, $3)); }
    | node
      { allocate($$, pt_parsetree_t(NODE_LIST_N, 1, NULL, $1)); }
    ;
node: name COLON distance
      { allocate($$, pt_parsetree_t(LEAF_N, 2, NULL, $1, $3)); }
    | LPAREN node_list RPAREN COLON distance
      { allocate($$, pt_parsetree_t(NODE_N, 2, NULL, $2, $5)); }
    ;
name: NAME { allocate($$, pt_parsetree_t(NAME_N, 0, strdup(yyget_text(context->scanner)))); }
    ;
distance:
      FLOAT
      {
        allocate($$, pt_parsetree_t(DISTANCE_N, 0, calloc(1, sizeof(double))));
        *((double *)$$->data) = atof(yyget_text(context->scanner));
      }
    ;
%%
